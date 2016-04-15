import pysam
import collections
import multiprocessing
from general_python import toolbox
from ngs_python.variant import annovar
import pandas as pd
from scipy.stats import fisher_exact

def pairDist(read1, read2):
    ''' Function calculates inner and outer distance for paired reads if
    they are determined to be concordant. Function takes two arguemnts:
    
    1) read1 - Pysam AlignedSegment for read1.
    2) read2 - Pysam AlignedSegment for read2.
    
    '''
    # Initialise return variable
    outer = inner = None
    # Process reads on same chromsome
    if read1.reference_id == read2.reference_id:
        # Process reads where read1 on reverse and read2 on forward strand
        if read1.is_reverse and not read2.is_reverse:
            # Check distance between potential pairs
            if (read1.reference_end >= read2.reference_end and
                read1.reference_start >= read2.reference_start):
                outer = read1.reference_end - read2.reference_start
                inner = read1.reference_start - read2.reference_end
        # Process reads where read1 on forward and read2 on reverse strand
        elif not read1.is_reverse and read2.is_reverse:
            # Check distance between potential pairs
            if (read2.reference_end >= read1.reference_end and
                read2.reference_start >= read1.reference_start):
                outer = read2.reference_end - read1.reference_start
                inner = read2.reference_start - read1.reference_end
    # Return return variable
    return(outer, inner)

def createChrDict(bamFile):
    ''' Function returns a dictionary of chromosome data where the
    key is the target ID (tid) of the chromosome and the values are
    a tuple of the chromosome name and the chromosome length. Function
    takes the following 1 argument:
    
    1)  bamFile - Open Pysam AlignmentFile
    
    '''
    # Create tuple of chrosomes
    chrList = zip(bamFile.references, bamFile.lengths)
    # Add to dictionary
    chrDict = {}
    for c in chrList:
        tid = bamFile.gettid(c[0])
        chrDict[tid] = c
    # Return dictionary
    return(chrDict)
    
def startInterval(read, intervalSize, chrDict):
    ''' Function returns interval around the start of a read. If the
    read is aligned to forward strand then the start is defined as the
    end of the read aligned to the most 5' portion of the reference.
    If the read is aligned to the reverse strand then the start is defined
    as end of the read aligned to the most 3' portion of the reference.
    The function takes the following 3 arguments:
    
    1)  read - Pysam AlignedSegment
    2)  intervalSize - Size of the returned interval
    3)  chrDict - A dictionary where the keys are the target IDs (tid)
        of the chromosomes and the values is a tuple of the chromsome
        name and length.
    
    Function returns a tuple of the region in BED format; using zero
    based indexing and open ended interval and where name is read name
    and score is read mapping quality.
    
    '''
    # Extract flag
    flag = read.flag
    # Find strand of read and extract read start position
    if flag & 16:
        strand = '-'
        site = read.reference_end
    else:
        strand = '+'
        site = read.reference_start
    # Calculate start and end of interval
    halfSize = int(intervalSize / 2)
    start = site - halfSize
    end = site + halfSize
    # Extract chromosome data and adjust start and end
    chrom, length = chrDict[read.reference_id]
    start = max(0, start)
    end = min(length, end)
    # Extract mapping quality
    mapq = read.mapping_quality
    # Extract read name
    if flag & 64:
        name = read.query_name + '/1'
    elif flag & 128:
        name = read.query_name + '/2'
    else:
        name = read.query_name
    # Return data
    return((chrom, start, end, name, mapq, strand))

def pairGenerator(bamFile, mapQ = 20):
    ''' Function extracts read pairs from  name sorted BAM files.
    Secondary alignments are ignored. Function takes two arguments:
    
    1)  bamFile - Name of input bamFile.
    2)  mapQ - Minimum mapping quality of reads.
    
    '''
    # Open bamFile and loop through
    readList = []
    currentName = ''
    inBam = pysam.AlignmentFile(bamFile, 'rb')
    while True:
        # Extract read data
        try:
            read = inBam.next()
            readName = read.query_name
        except StopIteration:
            readName = 'EndOfFile'
        # Process completed families
        if readName != currentName:
            # Count and process read-pairs
            if len(readList) == 2:
                read1, read2 = readList
                if read1.is_read1 and read2.is_read2:
                    yield(read1, read2)
                elif read1.is_read2 and read2.is_read1:
                    yield(read2, read1)
            # Reset read list and current name
            currentName = readName
            readList = []
        # Break loop at end of BAM file
        if readName == 'EndOfFile':
            break
        # Skip secondary alignments
        elif (read.flag & 256):
            continue
        # Skip reads below supplied mapping quality
        if read.mapping_quality <= mapQ:
            continue
        # Else append read to read list
        else:
            readList.append(read)

def baseCalls(
        read, groupdel = False
    ):
    # Extract cigar tuple and check
    cigartuple = read.cigartuples
    if cigartuple is None:
        return({})
    # Remove clipping from cigar
    cigartuple = [x for x in cigartuple if x[0] < 4]
    # Extract operations from tuple and check
    tupleSet = set([x[0] for x in cigartuple])
    if tupleSet.intersection([3, 6, 7, 8]):
        raise IOError('Irregular cigar found: %s' %(cigartuple))
    # Extract additional data
    sequence = list(read.query_alignment_sequence)
    quality = list(read.query_alignment_qualities)
    positions = read.get_reference_positions()
    # Process indels
    if tupleSet.intersection([1, 2]):
        # Trim clipping and indels from start of cigartuple
        while cigartuple:
            # Remove indels not flanked by mapped sequence
            if cigartuple[0][0] in [1,2]:
                sequence = sequence[cigartuple[0][1]:]
                quality = quality[cigartuple[0][1]:]
                cigartuple = cigartuple[1:]
            # Stop trimming
            else:
                break
        # Remove clipping and indels from end of read
        while cigartuple:
            # Remove indels not flanked by mapped sequence
            if cigartuple[-1][0] in [1,2]:
                sequence = sequence[:-cigartuple[-1][1]]
                quality = quality[:-cigartuple[-1][1]]
                cigartuple = cigartuple[:-1]
            # Stop trimming
            else:
                break
        # Initialise variables to extract indel data
        location = 0
        insertion = []
        deletion = []
        # Find position and length of indels
        if cigartuple:
            for element in cigartuple:
                # Process mapped elements of cigar
                if element[0] == 0:
                    location += element[1]
                # Process insertion elements of cigar
                elif element[0] == 1:
                    insertion.append((location - 1, location + element[1] ))
                # Process deletion elements of cigar
                elif element[0] == 2:
                    deletion.append((location, element[1] ))
        # Process insertions
        for i in insertion:
            # Concatenate insertions into single list element
            sequence[i[0]:i[1]] = ["".join(sequence[i[0]:i[1]])]
            # Alter quality to mean of bases
            quality[i[0]:i[1]] = [ sum(quality[i[0]:i[1]]) / 
                len(quality[i[0]:i[1]]) ]
        # Reverse and loop through deletion
        deletion.reverse()
        for d in deletion:
            # Group deletions
            if groupdel:
                # Alter sequence to include "-" to signify deletion
                sequence = sequence[:d[0]] + ['-' * d[1]] + sequence[d[0]:]
                # Alter quality to mean of flank
                dQual = (quality[d[0] - 1] + quality[d[0]]) / 2
                quality = quality[:d[0]] + [dQual] + quality[d[0]:]
                # Alter positions
                dPos = positions[d[0] - 1] + 1
                positions = positions[:d[0]] + [dPos] + positions[d[0]:]
            # Associate deletions with individual bases
            else:
                # Alter sequence to include "-" to signify deletion
                sequence = sequence[:d[0]] + ['-'] * d[1] + sequence[d[0]:]
                # Alter quality to mean of flank
                dQual = (quality[d[0] - 1] + quality[d[0]]) / 2
                quality = quality[:d[0]] + [dQual] * d[1] + quality[d[0]:]
                # Alter positions
                dPos = range(positions[d[0] - 1] + 1,
                    positions[d[0] - 1] + d[1] + 1)
                positions = positions[:d[0]] + dPos + positions[d[0]:]
    # Create and return output
    if len(sequence) == len(quality) and len(quality) == len(positions):
        positionDict = dict(zip(positions, zip(sequence, quality)))
        return(positionDict)
    else:
        raise IOError('Could not ascribe sequence to positions')

def extractPosition(
       openBam, chrom, position, minMapQ = 20, minBaseQ = 20, groupdel = False
    ):
    # check arguments
    toolbox.checkArg(chrom, 'str')
    toolbox.checkArg(position, 'int', mn = 1)
    toolbox.checkArg(minMapQ, 'int', mn = 0)
    toolbox.checkArg(minBaseQ, 'int', mn = 0)
    toolbox.checkArg(groupdel, 'bool')
    # Set variables for mapping
    mapQuality = []
    baseCounts = {}
    # Loop through reads covering BAM
    for read in openBam.fetch(chrom, position-1, position):
        # Store mapping quality and skip reads with low values
        mapQuality.append(read.mapping_quality)
        if mapQuality[-1] < minMapQ:
            continue
        # Extract base calls for read
        baseDict = baseCalls(read, groupdel)
        # Extract reads for base of interest or skip
        try:
            base, quality = baseDict[position - 1]
        except KeyError:
            continue
        # Skip bases of poor quality
        if quality < minBaseQ:
            continue
        # Add base to base dictionary
        if base in baseCounts:
            baseCounts[base][0] += 1
        else:
            baseCounts[base] = [1, 0]
        # Add forward strand count to base dictionary
        if not read.is_reverse:
            baseCounts[base][1] += 1
    # Calculate and return results
    if len(mapQuality) > 0:
        meanMap = sum(mapQuality) / len(mapQuality)
    else:
        meanMap = 0
    return(baseCounts, meanMap)

def extractVariantCountsProcess(
        variantList, bamFile, pipe, minMapQ = 20, minBaseQ = 20,
        groupdel = False
    ):
    ''' Function calculates the frequency at which specified nucleotides
    are found at specific chromosomal regions. Function takes 6 arguments:
    
    1)  variantList - a tuple/list of tuples/lists that contain four
        elements; the chromosome, position, reference, and variant e.g
        [('chr1', 1, 'A', 'T'), ('chr2', 100, 'C', 'G')].
    2)  bamFile - Full path to BAM file
    3)  minMapQ - Minimum mapping quality of read to extract base.
    4)  MinBaseQ - Minimum base quality to extract base.
    
    Function returns a pandas dataframe containing the following five columns:
    
    1)  varcount - Count of the variant calls.
    2)  refcount - Count of the reference calls.
    3)  varfor - Frequency of reference reads on forward strand.
    4)  reffor - Frequency of variant reads on forward strand.
    5)  mapqual - Mean mapping score of ALL reads spanning the postion.
    
    '''
    # check arguments
    toolbox.checkArg(minMapQ, 'int', mn = 0)
    toolbox.checkArg(minBaseQ, 'int', mn = 0)
    toolbox.checkArg(groupdel, 'bool')
    # Create output dataframe
    variantNames = [':'.join(map(str, x)) for x in variantList]
    outData = pd.DataFrame(
        columns = ['allcount', 'allfor', 'refcount', 'reffor', 'varcount',
            'varfor', 'mapqual'],
        index = variantNames
    )
    # Open bamFile
    bam = pysam.AlignmentFile(bamFile)
    # Loop through variants and extract counts
    for name, (chrom, position, reference, variant) in zip(variantNames,
        variantList):
        # Extract base data
        baseCounts, mapqual = extractPosition(
            openBam = bam, chrom = chrom, position = position,
            minMapQ = minMapQ, minBaseQ = minBaseQ, groupdel = groupdel
        )
        # Count total bases
        allcount = sum([x[0] for x in baseCounts.values()])
        allfor = sum([x[1] for x in baseCounts.values()])
        # Extract reference counts
        if reference in baseCounts:
            refcount, reffor = baseCounts[reference]
        else:
            refcount = 0
            reffor = 0
        # Extract variant frequency
        if variant in baseCounts:
            varcount, varfor = baseCounts[variant]
        else:
            varcount = 0
            varfor = 0
        # Add data to output list
        outData.loc[name] = [allcount, allfor, refcount, reffor, varcount,
            varfor, mapqual]
    # Close bam
    bam.close()
    # Send data down pipe and close
    pipe.send(outData)
    pipe.close()

def calculateVariantMetrics(
        variantList, bamList, sampleNames, annovarPath, buildver, database,
        tempprefix, minMapQ = 20, minBaseQ = 20, groupdel = False 
    ):
    ''' Function calculates the frequency at which specified nucleotides
    are found at specific chromosomal regions. Function takes 6 arguments:
    
    1)  variantList - a tuple/list of tuples/lists that contain three
        elements. These three elements are the chromosome, position and
        variant of interest e.g [('chr1', 1, 'A'), ('chr2', 100, 'C')].
    2)  bamList - Full path to BAM files. It is presumed that the first
        file in the list is the reference BAM file.
    3)  sampleNames - A list/tuple of samples names. Must be one for every
        BAM files
    4)  minMapQ - Minimum mapping quality of read to extract base.
    5)  MinBaseQ - Minimum base quality to extract base.
    6)  groupdel - Boolean, whether to group deletions to 5' base.
    
    Function returns a list of tuples where each tuple contains the following
    four elements:
    
    1)  Count of the variant at the described position.
    2)  Count of the variant on the forward strand.
    3)  Frequency of base forward strand reads.
    4)  Mean mapping score of ALL reads spanning the described postion.
        This value is calculated from reads with mapping qualities and base
        qualities below the desired threshold.
    
    '''
    # Check arguments
    toolbox.checkArg(minMapQ, 'int', mn = 0)
    toolbox.checkArg(minBaseQ, 'int', mn = 0)
    toolbox.checkArg(groupdel, 'bool')
    # Check supplied Names
    if not isinstance(sampleNames, (list, tuple)):
        raise IOError('sampleNasmes must be a list or a tuple')
    if len(sampleNames) != len(bamList):
        raise IOError('Must be a sample name for each BAM')
    # Create output dataframe
    varnames = [':'.join(map(str,x)) for x in variantList]
    outputData = pd.DataFrame(index = varnames, columns = ['chr','pos','ref',
        'var','minp'])
    for x in varnames:
        chrom, position, ref, var = x.split(':')
        outputData.loc[x] = [chrom, int(position), ref, var, 1]
    # Create variables to store pipe and process data
    pipeList = []
    processList = []
    # Create process for each BAM
    for bamFile in bamList:
        # Create pipes
        pipeRecv, pipeSend = multiprocessing.Pipe(False)
        pipeList.append(pipeRecv)
        # Create and start process
        proc = multiprocessing.Process(
            target = extractVariantCountsProcess,
            args = (variantList, bamFile, pipeSend, minMapQ, minBaseQ,
                groupdel)
        )
        proc.start()
        processList.append(proc)
        # Close send pipe
        pipeSend.close()
    # Extract variants for each BAM
    for number, (name, pipe, process) in enumerate(
            zip(sampleNames, pipeList, processList)
        ):
        # Extract data from pipe and close pipe and process
        sampleData = pipe.recv()
        pipe.close()
        process.join()
        # select desired data
        selectSample = pd.DataFrame(index = varnames)
        selectSample['cov'] = sampleData['allcount']
        selectSample['freq'] = sampleData['varcount'] / sampleData['allcount']
        selectSample['mapq'] = sampleData['mapqual']
        # Calculate pvalue
        if number:
            pvalue = []
            for x in zip(
                zip(refVar, refNonVar), 
                zip(
                    sampleData['varcount'],
                    sampleData['allcount'] - sampleData['varcount']
                )
            ):
                pvalue.append(fisher_exact(x)[1])
            selectSample['pvalue'] = pvalue
        # Extract count for reference
        else:
            refVar = sampleData['varcount']
            refNonVar = sampleData['allcount'] - sampleData['varcount']
        # Rename tables columns and append to output
        selectSample.columns = [name + '_' + x for x in selectSample.columns]
        outputData = pd.concat([outputData, selectSample], axis = 1)
    # Calculate minium pvalue and sort
    pvalueIndex = [x.endswith('pvalue') for x in outputData.columns]
    minp = outputData.loc[:,pvalueIndex].min(1)
    outputData['minp'] = minp
    # Add annovar annotation and concat to dataframe
    geneAnno = annovar.geneAnno2DF(variantList = variantList,
        path = annovarPath, buildver = buildver, database = database,
        tempprefix = tempprefix)
    outputData = pd.concat([outputData, geneAnno], axis = 1)
    outputData.sort_values('minp', inplace = True)
    # Return data
    return(outputData)

#data = calculateVariantMetrics(
#    variantList = [('13',93096846,'A','C'),('7',59673519,'T','A')],
#    bamList = [
#        '/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/328909-NORM_dedup_realign_recal.bam',
#        '/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/328909-T1R2_dedup_realign_recal.bam'
#    ],
#    sampleNames = ['NORM', 'T1R2'],
#    annovarPath = '/farm/babs/redhat6/software/annovar_2015Jun17/annotate_variation.pl',
#    buildver = 'mm10',
#    database = '/farm/babs/redhat6/software/annovar_2015Jun17/mousedb/',
#    tempprefix = '/farm/scratch/rs-bio-lif/rabino01/Elza/annotemp'
#)
#print data
