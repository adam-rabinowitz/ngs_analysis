import collections
import multiprocessing
import pysam
import os
import pandas as pd

def bedgraph_coverage(intervals):
    ''' Function creates coverage from a supplied interval list.
    Function takes 1 argument:
    
    1)  intervals - Either a pair of 'multiprocessing.connection' objects
    or an iterable. The pipe or iterable should deliver 
    
    Function return a list of 3-element tuples where the first element is
    the start of the region, the second element is the end of the region
    and the third element is the coverage
    '''
    # Create iterator from pipes or supplied iterator
    try:
        parentConn, childConn = intervals
        parentConn.close()
        interval_iterator = iter(childConn.recv, None)
        return = childConn.send
    except:
        interval_iterator = intervals
    # Create change dictionary and add elements from iterable
    changeDict = collections.OrderedDict()
    for chrom, start, end in interval_iterator:
        if chrom in changeDict:
            changeDict[chrom][start] += 1
            changeDict[chrom][end] -= 1
        else:
            changeDict[chrom] = collections.defaultdict(int)
            changeDict[chrom][start] += 1
            changeDict[chrom][end] -= 1
    # Create variables to generate bedgraph data
    bedgraph = pd.DataFrame(columns = ['chrom', 'start', 'end', 'cov'])
    chromosomes = changeDict.keys()
    for chrom in chromosomes:
        # Create a series from dictionary and calculate coverage
        changeSeries = pd.Series(changeDict.pop(chrom))
        changeSeries = changeSeries[changeSeries != 0]
        coverage = changeSeries.cumsum()
        if coverage.iloc[-1] != 0:
            raise ValueError('Error in coverage calculation')
        # Create data frame and add to output
        changeDataFrame = pd.DataFrame({'chrom' : chrom,
            'start' : changeSeries.index[:-1], 'end' : changeSeries.index[1:],
            'cov' : coverage.iloc[:-1]}, columns = ['chrom', 'start', 'end',
            'cov'])
        bedgraph = pd.concat((bedgraph, changeDataFrame))
    # Return bedgraph
    return(bedgraph)

x = bedgraph_coverage([('x',15,25),('y',30,40),('x',10,20)])
print x

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

