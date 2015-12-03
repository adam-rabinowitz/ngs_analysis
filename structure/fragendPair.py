# Import required modules
from Bio import SeqIO
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import bisect
from general_functions import iohandle
import collections

def findFragendSites(fasta, resite):
    ''' Function creates FragendDict object. The object contains
    the location of all fragends for eachh strand of all
    chromosomes within a FASTA file.
    '''
    # Process restriction enzyme size and create output dictionary
    resite = resite.upper()
    frags = {'chr': [], 'rsesite': resite}
    # Create sequence object for resite and reverse complent
    standard = Seq(self.resite)
    revcomp = standard.reve
    # Open and parse fasta file
    fastaHandle = open(fasta)
    fastaData = SeqIO.parse(fastaHandle,'fasta')
    # Loop through fasta file and extract fragend information for each chromosome
    for fasta in fastaData:
        # Extract name and sequence
        fName, fSequence = str(fasta.id), str(fasta.seq).upper()
        frags['chr'].append(fName)
        # Add re sites to dictionary using 1 based index
        forward = nt_search(fSequence, standard)[1:]
        if forward:
            self.frags[(fName,'+')] = [x + len(resite) for x in forward]
        else:
            self.frags[(fName,'+')] = []
        reverse = nt_search(fSequence, revcomp)[1:]
        if reverse:
            self.frags[(fName,'-')] = [x + 1 for x in reverse]
        else:
            self.frags[(fName,'-')] = []
    # Close input file and return data
    fastaHandle.close()
    return(frags)

def downstream(readIn, fragDict, fragOut = []):
    ''' Function to find the associated fragend for a read. For
    reads aligned to the '+' strand the fragend will be downstream
    of the start of the read. For reads aligned to the '-' strand
    the fragend will be upstream of the start of the read. Function
    takes three arguments:

    1)  readIn - A file/pipe/list providing a list/tuple with the
    following information for an aligned read: chrom, start, end
    and strand. Start is the most 5' base in the genome to which
    the read is aligned. End is the most 3' base in the genome to
    which the read is aligned

    2)  fragDict - A fragend dictionary generated using the
        findFragendSites function.

    3) fragOut - A file/pipe/list to ehich the output data will be
    sent

    Function returns a series of tuples containing the following
    four pices of information:

    1)  chrom - Chromsome to which the read is aligned.

    2)  fragLoc - Location of the centre of the restriction enzyme
    site from which the associated fragend is derived.

    3)  strand - Strand of the fragend.

    4)  distance - Distance between the start of the read and the
    fragends centre. For a '+' strand read distance is measured
    from 5' end of read. For a '-' strand read distance is measured
    from the 3' end of the read.
    '''
    # Create input and output objects
    inData = iohandle.handlein(readIn)
    outData = iohandle.handleout(fragOut)
    # Extract resite
    resite = fragDict['resite']
    # Sequentially process in data
    while True:
        # Get additional input or break loop
        try:
            chrom, start, end, strand = inData.next()
        except EOFError:
            break
        # Process variables
        start = int(start)
        end = int(end)
        fragLoc = None
        distance = None
        # Check relative location of start and end
        if end < start:
            raise ValueError("End must be more 3' than the start")
        # Extract fragends for chromosome and strand
        try:
            fragends = self.frags[(chrom, strand)]
        except KeyError:
            raise IOError('No data for %s strand on %s chromosome' %(
                chrom, strand))
        # Find location of fragend for forward strand read
        if strand == '+':
            fragIndex = bisect.bisect_left(
                fragends,
                end
            )
            # Return default data for invalid index
            if fragIndex < len(fragends):
                # Find fragend location
                fragLoc = fragends[fragIndex] - (len(resite) / 2)
                # Find distance between read and fragend
                distance = (fragLoc - start) + 1
        # Find location of fragend for reverse strand read
        elif strand == '-':
            # Find potential index of downstream fragend
            fragIndex = bisect.bisect_right(
                fragends,
                start
            ) - 1
            # Return default data for invalid index
            if fragIndex >= 0:
                # Find fragend location
                fragLoc = fragends[fragIndex] + (len(resite) / 2)
                # Find distance between read and fragend
                distance = (end - fragLoc) + 1
        # Add output data
        outData.add((chrom, fragLoc, strand, distance))
    # Close IO and return data
    inData.close()
    return(outData.close())

def fragend_pairs(pairIn, fasta, resite ,maxDistance, pairOut = []):
    ''' Function identifies and reports upstream fragends for HiC read pairs. The
    function takes 5 arguments:
    
    1)  pairIn - Read apit input object
    2)  fasta - Genome fasta file
    3)  reSite - Restriction enzyme recognition sequence
    4)  maxDistance - Maximum acceptable distance between start of read and RE site.
    5)  pairOut - Name of output gzipped file containing fragend ligations.
    
    '''
    # Create fragend dictionary
    fragDict = findFragendSites(fasta, resite)
    # Create variable to find pairs
    pairCounts = collections.OrderedDict([
        ('total', 0),
        ('no fragend', 0),
        ('fragend too distant', 0),
        ('intrachromosomal', 0),
        ('interchromosomal', 0),
        ('fragend distance', []),
        ('ligation distance', [])
    ])
    # Generate input and output object
    inData = iohandle.handlein(pairIn)
    outData = iohandle.handleout(pairOut)
    # Loop through pairs
    while True:
        try:
            pairData = inData.next()
        except EOFError:
            break
        # Count entries
        pairCounts['total'] += 1
        # Create output containg fragend data
        output = downstream([pairData[0:4],pairData[4:8]], fragDict)
        # Skip reads without identified fragends
        if output[0][1] == None or output[1][1] == None:
            pairCounts['no fragend'] += 1
            continue
        # Add fragend distance data for pairs with fragends
        pairCounts['fragend distance'].extend([
            output[0][3],
            output[1][3] 
        ])
        # Count and skip reads too distant from the fragend
        if output[0][3] > maxDistance or output[1][3] > maxDistance:
            pairCounts['fragend too distant'] += 1
            continue
        # Save to file accepted ligation pairs
        outData.add(output[0][0:3] + output[1][0:3])
        # Count interchromosomal ligations 
        if output[0][0] != output[1][0]:
            pairCounts['interchromosomal'] += 1
            # Count intrachromosomal ligations and store distance
        else:
            pairCounts['intrachromosomal'] += 1
            pairCounts['ligation distance'].append(
                abs(output[0][1] - output[1][1])
            )
    # Close files and return data
    inData.close()
    return(outData.close(),pairCounts)
