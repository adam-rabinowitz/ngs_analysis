import pysam
import numpy

def concordant(read1, read2, maxSize = 2000):
    ''' Function takes two pysam AlignedSegment reads and determines
    if they are concordant. If the reads are concoradant then outer and
    inner distance btween the reads will be returned. Function takes 
    three arguemnts:

    1)  read1 - AlignedSegment of read1.
    2)  read2 - AlignedSegment of read2.
    3)  maxSize - maximum outer distance of concordant reads.

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
            if (read2.reference_end >= read.reference_end and
                read2.reference_start >= read1.reference_start):
                outer = read2.reference_end - read1.reference_start
                inner = read2_reference_start - read1.reference_end
    # Retrun return variable
    return(outer, inner)

def pairGeneratorPaired(bamFile, mapQ = ''):
    ''' Function extracts read pairs from paired-end aligned, name
    sorted BAM files. Secondary alignments are ignored. Function
    takes two arguments:

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
                print 'Pair!'
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
        elif mapQ and read.mapping <= mapQ:
            continue
        # Else append read to read list
        else:
            readList.append(read)

def tophatInsertMetrics(bamFile, maxSize = 2000):
    innerList = []
    for read1, read2 in pairGeneratorPaired(bamFile):
        outer, inner = concordant(
            read1 = read1,
            read2 = read2,
            maxSize = maxSize
        )
        if outer and outer <= maxSize:
            innerList.append(inner)
    # Calculate and return metrics
    length = len(innerList)
    mean = numpy.mean(innerList, dtype=int)
    stdv = numpy.std(innerList, dtype=int, ddof=1)
    return(length, mean, stdv)
            

