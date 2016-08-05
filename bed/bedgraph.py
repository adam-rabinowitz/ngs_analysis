import collections
import multiprocessing
import pysam
import os
import pandas as pd

def bedgraph_coverage(intervals, zero=False, lengths=None):
    ''' Function creates coverage from a supplied interval list.
    Function takes 1 argument:
    
    1)  intervals - Either a pair of 'multiprocessing.connection' objects
    or an iterable.

    2) zero
    
    Function return a list of 3-element tuples where the first element is
    the start of the region, the second element is the end of the region
    and the third element is the coverage
    '''
    # Check arguments
    if zero and lengths is None:
        raise ValueError('Chromosome lengths required for zero coverage')
    # Create iterator from pipes or supplied iterator
    try:
        parentConn, childConn = intervals
        parentConn.close()
        interval_iterator = iter(childConn.recv, None)
        pipe = True
    except Exception:
        interval_iterator = intervals
        pipe = False
    # Create dictionary to store coverage changes
    coverageChange = collections.OrderedDict()
    if lengths:
        for chrom in lengths.keys():
            coverageChange[chrom] = collections.defaultdict(int)
        for chrom, start, end in interval_iterator:
            coverageChange[chrom][start] += 1
            coverageChange[chrom][end] -= 1
    else:
        for chrom, start, end in interval_iterator:
            try:
                coverageChange[chrom][start] += 1
                coverageChange[chrom][end] -= 1
            except KeyError:
                coverageChange[chrom] = collections.defaultdict(int)
                coverageChange[chrom][start] += 1
                coverageChange[chrom][end] -= 1
    # Create variables to generate bedgraph data
    bedgraph = pd.DataFrame(columns = ['chrom', 'start', 'end', 'cov'])
    for chrom in coverageChange.keys():
        # Extract chromosome dictionary and remove zero change
        chromChange = coverageChange.pop(chrom)
        chromChange = {k:v for k, v in chromChange.items() if v}
        if lengths:
            chromLength = lengths[chrom]
            chromChange[0] = chromChange.get(0, 0)
            chromChange[chromLength] = chromChange.get(chromLength, 0)
        # Create series, check lengths and create dataframe
        changeSeries = pd.Series(chromChange)
        if lengths and changeSeries.iloc[-1] > lengths[chrom]:
            raise ValueError('Intervals extend beyong chromosome')
        changeData = pd.DataFrame(
            { 'chrom' : chrom,
              'start' : changeSeries.index[:-1],
              'end' : changeSeries.index[1:],
              'cov' : changeSeries.cumsum().iloc[:-1] },
            columns = ['chrom', 'start', 'end', 'cov'])
        # Remove zero coverage and add to output
        if not zero:
            changeData = changeData[changeData['cov'] != 0]
        bedgraph = pd.concat((bedgraph, changeData))
    # Return data
    bedgraph.index = range(len(bedgraph.index))
    if pipe:
        childConn.send(bedgraph)
        childConn.close()
    else:
        return(bedgraph)


def bedgraph_bam(
        bam, intervals=None, zeros = False, min_map=0, duplicates=True
    ):
    # Create ordered dictionary of chromosome lengths
    bam = pysam.AlignmentFile(bam)
    lengths = collections.OrderedDict(zip(bam.references, bam.lengths))
    names = {bam.gettid(v):v for v in bam.references}
    # Check all intervals
    if intervals is None:
        checkedIntervals = [(None, None, None)]
    else:
        checkedIntervals = []
        for chrom, start, end in intervals:
            chromLength = lengths[chrom]
            if start is None:
                start = 0
            if end is None:
                end = chromLength
            if (start < 0
                or end > chromLength
                or start > end
                or not isinstance(start, int)
                or not isinstance(end, int)):
                raise ValueError('Invalid interval {}:{}-{}'.format(
                    chrom, start, end))
        # Add interval to checked list
        checkedIntervals.append((chrom, start, end))
    # Create flag for removing reads
    if duplicates:
        flag = 4
    else:
        flag = 1028
    # Create process for building bed graph
    connections = multiprocessing.Pipe()
    process = multiprocessing.Process(
        target = bedgraph_coverage,
        args = (connections, zeros, lengths)
    )
    process.start()
    parentConn, childConn = connections
    childConn.close()
    # Loop through bam files
    for chrom, start, end in checkedIntervals:
        for read in bam.fetch(chrom, start, end):
            # Skip unwanted and poorly mapped reads
            if read.flag & flag:
                continue
            if read.mapping_quality < min_map:
                continue
            # Send data
            parentConn.send((
                names[read.reference_id],
                max(start, read.reference_start),
                min(end, read.reference_end)
            ))
    # Retreive and return data and close process
    parentConn.send(None)
    bedgraph = parentConn.recv()
    parentConn.close()
    return(bedgraph)

bg = bedgraph_bam(
    bam = '/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/310123-NORM_dedup_realign_recal.bam'
)
print bg

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

