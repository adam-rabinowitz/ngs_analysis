import pysam
import numpy as np
import pandas as pd
import multiprocessing
import collections
import datetime

class coverage(object):

    def __init__(self, bam):
        self.bam = bam
        bamFile = pysam.AlignmentFile(bam)
        self.length = collections.OrderedDict(
            zip(bamFile.references, bamFile.lengths))
        self.name = {bamFile.gettid(x):x for x in bamFile.references}
    
    def extract_intervals(
            self, intervals, pipe_out, map_quality = 0, remove_dup = False
        ):
        bamFile = pysam.AlignmentFile(self.bam)
        for interval in intervals:
            pipe_out.send(interval)
            chrom, start, end = interval
            for read in bamFile.fetch(chrom, start, end):
                if read.is_unmapped:
                    continue
                if read.mapping_quality < map_quality:
                    continue
                if remove_dup and read.is_duplicate:
                    continue
                pipe_out.send((
                    max(start, read.reference_start),
                    min(end, read.reference_end)
                ))
        pipe_out.close()
        bamFile.close()
    
    def calculate_coverage_change(
            self, pipe_in, pipe_out
        ):
        region = None
        while True:
            try:
                data = pipe_in.recv()
            except EOFError:
                if region is not None:
                    pipe_out.send((region, coverage_change))
                break
            try:
                start, stop = data
                coverage_change[start] += 1
                coverage_change[stop] -= 1
            except ValueError:
                if region is not None:
                    pipe_out.send((region, coverage_change))
                region = data
                coverage_change = collections.defaultdict(int)
        pipe_in.close()
        pipe_out.close()
    
    def generate_bedgraph(
            self, pipe_in, pipe_out, zero_cov = False
        ):
        while True:
            try:
                interval, coverage_change = pipe_in.recv()
            except EOFError:
                break
            chrom, start, end = interval
            coverage_change = {k:v for k, v in coverage_change.items() if v}
            coverage_change[start] = coverage_change.get(start, 0)
            coverage_change[end] = coverage_change.get(end, 0)
            change_series = pd.Series(coverage_change)
            change_df = pd.DataFrame(
                {'chrom' : chrom, 'start' : change_series.index[:-1],
                    'end' : change_series.index[1:],
                    'cov' : change_series.cumsum().iloc[:-1]},
                columns = ['chrom', 'start', 'end', 'cov']
            )
            if not zero_cov:
                change_df = change_df[change_df['cov'] != 0]
            change_df.index = range(len(change_df.index))
            pipe_out.send((interval, change_df))
        pipe_in.close()
        pipe_out.close()

    def generate_interval_bedgraph(
            self, intervals = None, map_quality = 0, remove_dup = False,
            zero_cov = False
        ):
        # Process intervals
        if intervals is None:
            intervals = []
            for chrom in self.length:
                intervals.append((chrom, 0, self.length[chrom]))
        start_time = datetime.datetime.now()
        # Start process extracting intervals
        conn1Recv, conn1Send = multiprocessing.Pipe(False)
        process1 = multiprocessing.Process(
            target = self.extract_intervals,
            args = (intervals, conn1Send, map_quality, remove_dup)
        )
        process1.start()
        conn1Send.close()
        # Start process building coverage
        conn2Recv, conn2Send = multiprocessing.Pipe(False)
        process2 = multiprocessing.Process(
            target = self.calculate_coverage_change,
            args = (conn1Recv, conn2Send)
        )
        process2.start()
        conn1Recv.close()
        conn2Send.close()
        # Start process generating bedgraph
        conn3Recv, conn3Send = multiprocessing.Pipe(False)
        process3 = multiprocessing.Process(
            target = self.generate_bedgraph,
            args = (conn2Recv, conn3Send, zero_cov)
        )
        process3.start()
        conn2Recv.close()
        conn3Send.close()
        # Extract data
        output = pd.DataFrame(columns = ['chrom', 'start', 'end', 'cov'])
        while True:
            try:
                interval, bedgraph = conn3Recv.recv()
            except EOFError:
                break
            print interval
            print bedgraph
            print str(datetime.datetime.now() - start_time)
        conn3Recv.close()
        for process in [process1, process2, process3]:
            process.join()
    
    
            

x = coverage('/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/310123-NORM_dedup_realign_recal.bam')
x.generate_interval_bedgraph()



#def create_bins(sizes, width, fixed = True, overlap = 0, index = 0):
#    # Check index argument
#    if not index in (0, 1):
#        raise ValueError('index argument must be 0 or 1')
#    # Create output object and calculate span
#    bins = pd.DataFrame(columns = ['chrom','start','end'])
#    span = width - overlap
#    count = 0
#    # Loop through chromsome and generate bin data frame
#    for chrom, length in sizes:
#        length = float(length)
#        # Generate bin start and end for equal width bins
#        if fixed:
#            # Calculate end of first bin
#            remainder = (length - width) % span
#            end = np.floor(remainder/2) + width
#            # Calculate start and end of all bins
#            binEnd = np.arange(end, length + 1, span, dtype = np.uint32)
#            binStart = binEnd - width
#        # Generate bin start and end sites for unequal width bins
#        else:
#            # Calculate number of bins
#            binNo = np.ceil((length - width) / span) + 1
#            # Calculate excess number of elements within bins
#            excess = (width + (span * (binNo - 1))) - length
#            # Calculate bin size 
#            largeBin = int(width - np.floor(excess / binNo))
#            smallBin = int(width - np.ceil(excess / binNo))
#            # Calculate bin numbers
#            smallBinNo = int(excess % binNo)
#            largeBinNo = int(binNo - smallBinNo)
#            # Generate bin width
#            binWidth = np.array([largeBin] * largeBinNo + 
#                [smallBin] * smallBinNo, dtype = np.uint32)
#            # Generate bins
#            binEnd = np.cumsum(binWidth - overlap) + overlap
#            binStart = binEnd - binWidth
#        # Create chromosome dataframe
#        chromBins = pd.DataFrame(
#            {'chrom':chrom, 'start':binStart, 'end':binEnd},
#            index = np.arange(count, count + len(binStart))
#        )
#        if index == 1:
#            chromBins['start'] += 1
#        # Adjust bin count and store dataframe
#        count += len(binStart)
#        bins = pd.concat((bins, chromBins), axis=0)
#    # Return dataframe
#    bins = bins[['chrom','start','end']]
#    return(bins)
#
#def create_bins_bam(
#    bam, width, fixed = True, overlap = 0, index = 1
#):
#    # Create chromsome list
#    chrList = []
#    bam = pysam.AlignmentFile(bam)
#    for chrTuple in zip(bam.references, bam.lengths):
#        chrList.append(chrTuple)
#    bam.close()
#    # Create and return bed file
#    bins = create_bins(sizes = chrList, width = width, fixed = fixed,
#        overlap = overlap, index = index)
#    return(bins)
#    
#def count_reads_region(
#        bam, chrom, start, end, overlap = 'any', mapq = 0, index = 0
#    ):
#    # Check index and adjust start coordinates
#    if not index in (0, 1):
#        raise ValueError('Index must be 0 or 1')
#    start -= index
#    # Initialize counter
#    counter = 0
#    # Loop through counter
#    for read in bam.fetch(chrom, start, end):
#        # Skip unmapped and poorly mapped reads
#        if read.is_unmapped:
#            continue
#        if read.mapping_quality < mapq:
#            continue
#        # Count if overlap is 'any'
#        if overlap == 'any':
#            counter += 1
#        # Process alternative overlaps
#        else:
#            # Extract position
#            if overlap == 'start':
#                if read.is_reverse:
#                    position = read.reference_end
#                else:
#                    position = read.reference_start
#            elif overlap == 'end':
#                if read.is_reverse:
#                    position = read.reference_start
#                else:
#                    position = read.reference_end
#            else:
#                raise ValueError('Unrecognised overlap argument')
#            # Count if position is in desired interval
#            if start <= position <= end:
#                counter += 1
#    # Return counter
#    return(counter)
#
#def count_reads_bam(
#        bam, bed, overlap = 'any', index = 0, mapq = 0, pipe = None
#    ):
#    # check bed indices are unique
#    if len(bed.index) > len(set(bed.index)):
#        raise ValueError('bed indices must be unique')
#    # Open bam file and create index
#    bam = pysam.AlignmentFile(bam)
#    output = pd.Series(index = bed.index, dtype = np.int32)
#    # Loop through bed and extract read counts
#    for dfindex, chrom, start, end in bed.itertuples():
#        output[dfindex] = count_reads_region(bam = bam, chrom = chrom,
#            start = start, end = end, overlap = overlap, mapq = mapq, 
#            index = index)
#    # Close bam and return values
#    bam.close()
#    if pipe:
#        pipe[0].close()
#        pipe[1].send(output)
#        pipe[1].close()
#    else:
#        return(output)
#
#def count_reads_multibam(
#        bamList, bed, overlap = 'any', index = 0, mapq = 0
#    ):
#    # Create output dataframe
#    output = pd.DataFrame(index = bed.index, columns = bamList,
#        dtype = np.int32)
#    # Create and store process for each bam
#    processList = []
#    for bam in bamList:
#        pipe = multiprocessing.Pipe()
#        process = multiprocessing.Process(
#            target = count_reads_bam,
#            args = (bam, bed, overlap, index, mapq, pipe)
#        )
#        process.start()
#        pipe[1].close()
#        processList.append((bam, process, pipe[0]))
#    # Extract data for each bam
#    for bam, process, pipe in processList:
#        data = pipe.recv()
#        pipe.close()
#        process.join()
#        output[bam] = data
#    # Return output
#    return(output)
#
#bamList = [
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557175_sort.bam',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557176_sort.bam',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557178_sort.bam',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1559300_sort.bam',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1564297_sort.bam',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1564298_sort.bam',
#]
#
#bed = create_bins_bam(
#    bam = bamList[0],
#    width = 1000,
#    fixed = True,
#    index = 1
#)
#data = count_reads_multibam(
#    bamList = bamList,
#    bed = bed,
#    overlap = 'end',
#    mapq = 0,
#    index = 1
#)
#data = pd.concat((bed,data), axis = 1)
#data.to_csv(
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/binCount.txt',
#    sep = '\t',
#    header = True,
#    index = False,
#    float_format = '%.0f'
#)
