# Import required modules
from Bio import SeqIO
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import bisect
import collections
import gzip
from general_functions import writeFile

def findFragendSites(fasta, resite):
    ''' Function creates FragendDict object. The object contains
    the location of all fragends for eachh strand of all
    chromosomes within a FASTA file.
    '''
    # Process restriction enzyme size and create output dictionary
    resite = resite.upper()
    frags = {'resite': resite}
    # Create sequence object for resite and reverse complent
    standard = Seq(resite)
    revcomp = standard.reverse_complement()
    # Open and parse fasta file
    fastaHandle = open(fasta)
    fastaData = SeqIO.parse(fastaHandle,'fasta')
    # Loop through fasta file and extract fragend information for each chromosome
    for fasta in fastaData:
        # Extract name and sequence
        fName, fSequence = str(fasta.id), str(fasta.seq).upper()
        # Add re sites to dictionary using 1 based index
        forward = nt_search(fSequence, standard)[1:]
        if forward:
            frags[(fName,'+')] = [x + len(resite) for x in forward]
        else:
            frags[(fName,'+')] = []
        reverse = nt_search(fSequence, revcomp)[1:]
        if reverse:
            frags[(fName,'-')] = [x + 1 for x in reverse]
        else:
            frags[(fName,'-')] = []
    # Close input file and return data
    fastaHandle.close()
    return(frags)

def downstream(readIn, fragDict):
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
    # Extract resite and create output variable
    resite = fragDict['resite']
    output = []
    # Sequentially process in data
    for read in readIn:
        # Process variables
        chrom, start, end, strand = read
        start = int(start)
        end = int(end)
        fragLoc = None
        distance = None
        # Check relative location of start and end
        if end < start:
            raise IOError("End must be 'right' of the start")
        # Extract fragends for chromosome and strand
        try:
            fragends = fragDict[(chrom, strand)]
        except KeyError:
            raise IOError('No data for %s strand on %s chromosome' %(
                strand, chrom))
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
        if fragLoc == None:
            output.append(None)
        else:
            output.append((chrom, fragLoc, strand, distance))
    # Close IO and return data
    return(output)

def fragendPairs(pairIn, fasta, resite ,maxDistance, fragendOut):
    ''' Function identifies and reports upstream fragends for HiC read pairs. The
    function takes 5 arguments:
    
    1)  pairIn - Read apit input object
    2)  fasta - Genome fasta file
    3)  reSite - Restriction enzyme recognition sequence
    4)  maxDistance - Maximum acceptable distance between start of read and RE site.
    5)  pairOut - Name of output gzipped file containing fragend ligations.
    
    '''
    # Create fragend dictionary and metrics dictionary
    fragDict = findFragendSites(fasta, resite)
    fragendCounts = collections.defaultdict(int)
    fragendCounts['fragDist'] = []
    fragendCounts['ligDist'] = []
    # Open input and output file
    if pairIn.endswith('.gz'):
        inFile = gzip.open(pairIn, 'r')
    else:
        inFile = open(pairIn, 'r')
    outFile = writeFile.writeFileProcess(fragendOut)
    for pair in inFile:
        pair = pair.strip().split('\t')
        # Count entries
        fragendCounts['total'] += 1
        # Create output containg fragend data
        output = downstream([pair[0:4],pair[4:8]], fragDict)
        # Skip reads without identified fragends
        if output[0] == None or output[1] == None:
            fragendCounts['none'] += 1
            continue
        # Add fragend distance data for pairs with fragends
        fragendCounts['fragDist'].extend([
            output[0][3],
            output[1][3] 
        ])
        # Count and skip reads too distant from the fragend
        if output[0][3] > maxDistance or output[1][3] > maxDistance:
            fragendCounts['distant'] += 1
            continue
        # Save to file accepted ligation pairs
        outData = '\t'.join(map(str,output[0][0:3] + output[1][0:3]))
        outFile.add(outData + '\n')
        # Count interchromosomal ligations 
        if output[0][0] != output[1][0]:
            fragendCounts['interchromosomal'] += 1
            # Count intrachromosomal ligations and store distance
        else:
            fragendCounts['intrachromosomal'] += 1
            fragendCounts['ligDist'].append(
                abs(output[0][1] - output[1][1])
            )
    # Close files and return data
    inFile.close()
    outFile.close()
    return(fragendCounts)
