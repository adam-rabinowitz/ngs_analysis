import collections
import pysam
import os

def intlist2covlist(intlist):
    ''' Function creates coverage from a supplied interval list.
    Function takes 1 argument:
    
    1)  intlist - A list of 2-element tuples/lists where the first element
        is the start of the interval and the second element if the end of
        the interval.
    
    Function return a list of 3-element tuples where the first element is
    the start of the region, the second element is the end of the region
    and the third element is the coverage
    '''
    # Create change dictionary
    changeDict = collections.defaultdict(int)
    for start, end in intlist:
        changeDict[start] += 1
        changeDict[end] -= 1
    # Extract keys from dictionary and sort
    changePositions = changeDict.keys()
    changePositions.sort()
    # Set variables to calculate coverage
    coverage = 0
    coverageStart = 0
    coverageList = []
    # Loop through changes and calculate coverage
    for position in changePositions:
        # Extract change and skip zeros
        change = changeDict[position]
        if change == 0:
            continue
        # Write output if coverage != 0
        if coverage != 0:
            # Append region to list
            coverageList.append((coverageStart, position, coverage))
        # Adjust variables
        coverage += change
        coverageStart = position
    # Return coverage list
    return(coverageList)


def intlist2bg(pipe, bedgraph):
    ''' Function to parse list of overlapping genomic intervals and
    create a bedgraph file. Function takes two arguments.
    
    1)  pipe - Python pipe down which the chromsome name and interval list
        will be sent.
    2)  bedgraph - Name of bedgraph file.
    
    '''
    # Print pid
    print 'intList2bg', os.getpid()
    # Open out file
    with open(bedgraph, 'w') as bg:
        # Loop to process lists
        while True:
            # Extract data from pipe
            inData = pipe.recv()
            # Create change dictionary from input list
            if isinstance(inData, list):
                coverageList = intlist2covlist(inData)
            # Set chromosome name from input string
            elif isinstance(inData, str):
                chromosome = inData
                continue
            # Break loop for input None
            elif inData is None:
                pipe.close()
                break
            # Raise error for unexpected type
            else:
                print inData
                raise TypeError('Unexpected entry %s')
            # Print calculate coverage to file
            for start, end, coverage in coverageList:
                bg.write('%s\t%s\t%s\t%s\n' %(chromosome, start, end, coverage))

def bam2bg(bam, bedgraph):
    ''' Function extracts all mapped reads a BAM file and creates a bedgrpah
    file containing absolute coverage. Function takes two arguments:
    
    1)  bam - Full path to input bam file.
    2)  bedgraph - Full path to output bedgraph file.
    
    '''
    # Print pid
    print 'bam2bg', os.getpid()
    # Open bam file and extract chromosome names
    bamFile = pysam.AlignmentFile(bam)
    chromosomes = bamFile.references
    # Create pipe
    pipeRecv, pipeSend = multiprocessing.Pipe(False)
    # Create and start int2cov process
    intlist2bgProcess = multiprocessing.Process(
        target = intlist2bg,
        args = (pipeRecv, bedgraph)
    )
    intlist2bgProcess.start()
    # Close unwanted pipes
    pipeRecv.close()
    # Loop through chromsome names and extract intervals
    for chrom in chromosomes:
        # Send chromosome name
        pipeSend.send(chrom)
        # Create counter and iterator
        count = 0
        chromIterator = bamFile.fetch(chrom)
        # Create initial values from first mapped read
        while True:
            # Extract read
            try:
                read = chromIterator.next()
            except StopIteration:
                break
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Extract start and end, create interval list and maximum
            start, end = read.reference_start, read.reference_end
            intervalList = [(start, end)]
            maxList = end
            # Count read and break
            count += 1
            break
        # Loop through remaining reads in chromsome
        for read in chromIterator:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Count read and extract start and end
            count += 1
            start, end = read.reference_start, read.reference_end
            # Process non-overlapping interval
            if start > maxList:
                # Send current chromosome interval list
                pipeSend.send(intervalList)
                # Create new interval list and list maximum
                intervalList = [(start, end)]
                maxList = end
            # Process overlapping intervals
            else:
                # Add interval to list and recalculate list maximum
                intervalList.append((start, end))
                maxList = max(maxList, end)
        # Send last list and print chromosome metrics
        if len(intervalList) > 0:
            pipeSend.send(intervalList)
        print '%s: %s' %(chrom, count)
    # Add poisin pill, close pipes and join process
    pipeSend.send(None)
    pipeSend.close()
    intlist2bgProcess.join()

import multiprocessing
bam2bg('/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/310123-NORM_dedup_realign_recal.bam', '/farm/home/rabino01/test.bedgraph')
