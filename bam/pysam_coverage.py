import os
import pysam
import numpy as np
import pandas as pd
import multiprocessing
import collections

class single_coverage(object):
    
    def __init__(self, bam):
        self.bam = bam
        bamFile = pysam.AlignmentFile(bam)
        self.length = collections.OrderedDict(
            zip(bamFile.references, bamFile.lengths))
        bamFile.close()
    
    def __extract_intervals(
            self, intervals, conn_out, map_quality, remove_dup,
            remove_secondary
        ):
        # Process connections
        conn_out[0].close()
        # Set filter flag to skip unwanted reads
        filter_flag = 516
        if remove_dup:
            filter_flag += 1024
        if remove_secondary:
            filter_flag += 256
        # Open bam file and loop through intervals
        bamFile = pysam.AlignmentFile(self.bam)
        for interval in intervals:
            # Send and extract interval data
            conn_out[1].send(interval)
            chrom, start, end = interval
            # Extract reads for interval
            for read in bamFile.fetch(chrom, start, end):
                # Skip unmapped, failed, duplicate and secondary reads
                if read.flag & filter_flag:
                    continue
                # Skip reads below specified mapping quality
                if read.mapping_quality < map_quality:
                    continue
                # Extract blocks for read
                conn_out[1].send(read.get_blocks())
        # Tidy up
        bamFile.close()
        conn_out[1].close()
    
    def __coverage_change(
            self, intervals, conn_out, map_quality, remove_dup,
            remove_secondary
        ):
        # Process output connection and create input connection
        conn_out[0].close()
        conn_in = multiprocessing.Pipe(False)
        # Start process extracting intervals
        process = multiprocessing.Process(
            target = self.__extract_intervals,
            args = (intervals, conn_in, map_quality, remove_dup,
                remove_secondary)
        )
        process.start()
        conn_in[1].close()
        # Process initial interval
        cov_change = collections.defaultdict(int)
        interval = conn_in[0].recv()
        interval_start, interval_end = interval[1:]
        conn_out[1].send(interval)
        current_end = 0
        # Loop through remaining intervals
        while True:
            # Extract data or return final coverage dictionary and break
            try:
                data = conn_in[0].recv()
            except EOFError:
                conn_out[1].send(cov_change)
                break
            # Process read blocks
            if isinstance(data[0][0], int):
                # Send coverage change if new read is non-overlapping
                if current_end and data[0][0] > current_end:
                    conn_out[1].send(cov_change)
                    cov_change = collections.defaultdict(int)
                current_end = max(current_end, data[-1][-1])
                # Add data to coverage change object
                for block_start, block_end in data:
                    if block_end <= interval_start:
                        continue
                    if block_start >= interval_end:
                        continue
                    cov_change[max(interval_start, block_start)] += 1
                    cov_change[min(interval_end, block_end)] -= 1
            # Process new interval
            else:
                # Process coverage change object
                conn_out[1].send(cov_change)
                cov_change = collections.defaultdict(int)
                interval = data
                interval_start, interval_end = interval[1:]
                conn_out[1].send(interval)
                current_end = 0
        # Clean up processes and pipes
        conn_in[0].close()
        process.join()
        conn_out[1].close()
    
    def __coverage_array(
            self, intervals, conn_out, map_quality, remove_dup,
            remove_secondary
        ):
        # Proces output connection and create input connection
        conn_out[0].close()
        conn_in = multiprocessing.Pipe(False)
        # Start process extracting coverage change
        process = multiprocessing.Process(
            target = self.__coverage_change,
            args = (intervals, conn_in, map_quality, remove_dup,
                remove_secondary)
        )
        process.start()
        conn_in[1].close()
        # Extract data from pipe
        while True:
            # Extract data or break
            try:
                data = conn_in[0].recv()
            except EOFError:
                break
            # Calculate coverage intervals or break
            try:
                cov_change = {k:v for k,v in data.items() if v}
            except AttributeError:
                conn_out[1].send(data)
                continue
            # Convert to numpy array, sort, calculate coverage and send
            cov_array = np.array(cov_change.items(),
                dtype=[('pos', np.int32),( 'cov', np.int32)])
            if cov_array.size:
                cov_array = cov_array[cov_array['pos'].argsort()]
                cov_array['cov'] = np.cumsum(cov_array['cov'])
                conn_out[1].send(cov_array)
        # Clean up pipes and processes
        conn_in[0].close()
        process.join()
        conn_out[1].close()
    
    def __start_coverage_processes(
            self, intervals, map_quality, remove_dup, remove_secondary
        ):
        # Check intervals
        self.check_intervals(intervals)
        # Check other arguments
        if not isinstance(map_quality, int):
            raise TypeError('map_quality must be integer')
        if map_quality < 0:
            raise ValueError('map_quality must be non-negative')
        if not isinstance(remove_dup, bool):
            raise TypeError('remove_dup must be a bool')
        if not isinstance(remove_secondary, bool):
            raise TypeError('remove_secondary must be a bool')
        # Start process building coverage
        conn_in = multiprocessing.Pipe(False)
        process = multiprocessing.Process(
            target = self.__coverage_array,
            args = (intervals, conn_in, map_quality, remove_dup,
                remove_secondary)
        )
        process.start()
        conn_in[1].close()
        # Return process and pipe
        return((process, conn_in[0]))
    
    def check_intervals(self, intervals):
        # Check intervals is a list/tuple of length > 0
        if not isinstance(intervals, (list, tuple)):
            raise TypeError('intervals must be a list or tuple')
        if not len(intervals) > 0:
            raise ValueError('intervals must have length > 0')
        # Check all intervals
        for chrom, start, end in intervals:
            # Check chromsome
            if not isinstance(chrom, str):
                raise TypeError('Interval chromosome must be a string')
            try:
                chromLength = self.length[chrom]
            except KeyError:
                raise ValueError('Invalid chromosome: {}'.format(chrom))
            # Check start
            if not isinstance(start, int):
                raise TypeError('Interval starts must be integers')
            if not start >= 0:
                raise ValueError('Invalid start: {}'.format(start))
            # Check end
            if not isinstance(end, int):
                raise TypeError('Interval ends must be integers')
            if not end <= chromLength:
                raise ValueError('Invalid end: {}:{}'.format(chrom, end))
            # Check interval
            if start >= end:
                raise ValueError('Invalid interval: {} {}'.format(start, end))
    
    def mean_coverage_each(
            self, intervals, map_quality = 0, remove_dup = False,
            remove_secondary = False
        ):
        # Start processes to calculate coverage arrays
        process, pipe = self.__start_coverage_processes(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup, remove_secondary=remove_secondary)
        # Create output series
        outArray = np.zeros(len(intervals), dtype=np.float64)
        # Sequentially extract data from pipe
        count = 0
        data = pipe.recv()
        covDict = collections.defaultdict(int)
        covDict[0] = data[2] - data[1]
        while True:
            # Extract data from pipe or terminate loop
            try:
                data = pipe.recv()
            except EOFError:
                covArray = np.array(covDict.items(), dtype=np.int32)
                meanCov = np.average(covArray[:,0], weights=covArray[:,1])
                outArray[count] = meanCov
                break
            # Extract lengths
            try:
                lengths = np.diff(data['pos'])
            # Or process old interval
            except TypeError:
                # Calculate mean coverage
                covArray = np.array(covDict.items(), dtype=np.int32)
                meanCov = np.average(covArray[:,0], weights=covArray[:,1])
                outArray[count] = meanCov
                # Set parameters for new interval
                count += 1
                covDict = collections.defaultdict(int)
                covDict[0] = data[2] - data[1]
                continue
            # Add lengths and coverage to output
            for length, coverage in zip(lengths, data['cov']):
                covDict[coverage] += length
                covDict[0] -= length
        # Clean up processes and pipes
        pipe.close()
        process.join()
        return(outArray)
    
    def coverage_count_all(
            self, intervals, map_quality = 0, remove_dup = False,
            remove_secondary = False, max_cov = np.inf
    ):
        # Check argument
        if not np.isinf(max_cov) and not isinstance(max_cov, int):
            raise TypeError('max_cov must be infinite or integer')
        if not max_cov > 0:
            raise ValueError('max_cov > 0')
        # Start processes to calculate coverage arrays
        process, pipe = self.__start_coverage_processes(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup, remove_secondary=remove_secondary)
        # Extract data
        covCount = collections.defaultdict(np.int32)
        while True:
            # Extract data from pipe or terminate loop
            try:
                data = pipe.recv()
            except EOFError:
                break
            # Extract lengths
            try:
                lengths = np.diff(data['pos'])
            # Or process old interval
            except TypeError:
                covCount[0] += (data[2] - data[1])
                continue
            # Add lengths and coverage to output
            for length, coverage in zip(lengths, data['cov']):
                covCount[min(coverage, max_cov)] += length
                covCount[0] -= length
        # Clean up processes and pipes and return data
        pipe.close()
        process.join()
        return(covCount)
    
    def coverage_histogram(
            self, intervals, max_cov = 1000, map_quality = 0,
            remove_dup = False, remove_secondary = False
        ):
        # Extract coverage count dictionary
        covCount = self.coverage_count_all(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup, remove_secondary=remove_secondary,
            max_cov=max_cov
        )
        # Calculate histogram
        counts = np.zeros(max_cov + 1, dtype=np.int32)
        for c in xrange(max_cov + 1):
            counts[c] = covCount[c]
        revCumSum = counts[::-1].cumsum()
        totalBases = counts.sum(dtype=np.float64)
        covHist = revCumSum[::-1] / totalBases
        return(covHist)
    
    def mean_coverage_all(
            self, intervals, map_quality = 0, remove_dup = False,
            remove_secondary = False
        ):
        # Extract coverage count dictionary
        covCount = self.coverage_count_all(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup, remove_secondary=remove_secondary,
            max_cov=np.inf
        )
        # Calculate histogram
        covCounts = np.array(covDict.items(), dtype=np.int32)
        meanCov = np.average(covCounts)
        totalBases = covArray[:,1].sum(dtype=np.float64)
        covHist = revCumSum / totalBases
        return(covHist)
    
    def bedgraph_file(
            self, bedgraph, intervals, min_cov = 0, map_quality = 0,
            remove_dup = False, remove_secondary = False
        ):
        # Check argument
        if not isinstance(min_cov, int):
            raise TypeError('min_cov must be integer')
        if min_cov < 0:
            raise ValueError('min_cov cannot be negative')
        # Open output file
        outFile = open(bedgraph, 'w')
        # Start processes to calculate coverage arrays
        process, pipe = self.__start_coverage_processes(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup, remove_secondary=remove_secondary)
        # Extract data
        interval_end = current_end = 0
        while True:
            # Extract data from pipe or terminate loop
            try:
                data = pipe.recv()
            except EOFError:
                if interval_end > current_end and min_cov == 0:
                    outFile.write('{}\t{}\t{}\t0\n'.format(
                        chrom, current_end, interval_end))
                break
            # Generate an array of the coverage
            try:
                covArray = np.array(
                    zip(data['pos'][:-1], data['pos'][1:], data['cov'][:-1]),
                    dtype = [('start', np.int32), ('end', np.int32),
                        ('cov', np.int32)])
            # Or process old interval
            except TypeError:
                if interval_end > current_end and min_cov == 0:
                    outFile.write('{}\t{}\t{}\t0\n'.format(
                        chrom, current_end, interval_end))
                chrom, current_end, interval_end = data
                continue
            # Add zero-coverage gaps to output
            if covArray['start'][0] > current_end and min_cov == 0:
                outFile.write('{}\t{}\t{}\t0\n'.format(
                    chrom, current_end, covArray['start'][0]))
            # Add coverage to output
            for start, end, cov in covArray:
                if cov >= min_cov:
                    outFile.write('{}\t{}\t{}\t{}\n'.format(
                        chrom, start, end, cov))
            # Set current end
            current_end = covArray['end'].max()
        # Clean up processes and pipes and files
        outFile.close()
        pipe.close()
        process.join()
    
    # Function to create genome bins
    def create_bins(self, binSize, binEqual):
        ''' A function to create bins covering the genome.
        
        Args:
            binsize (int)- Maximum bin size.
            binequal (bool)- If True then all bins will be of an equal size.
                and if a chromosome is not a multiple of the binsize te most
                telomeric regions will be ignored. If False then bin sizes
                may be reduced to cover the chromosome as equally as possible.
        
        Returns:
            binarray -A numpy array of the chromosome, start site and end of
                the bins. Chromosome positions are 0-indexed half-open.
        
        '''
        # Check arguments
        if not isinstance(binSize, int):
            raise TypeError('binSize must be integer')
        if not binSize > 0:
            raise ValueError('binSize must be >0')
        if not isinstance(binEqual, bool):
            raise TypeError('binEqual must be bool')
        # Create output data variables
        binDict = collections.OrderedDict()
        count = 0
        # Loop through chromsome and generate bin data frame
        for chrName, chrLength in self.length.items():
            # Generate bin start and end for equal width bins
            if binEqual:
                # Calculate initial start and end
                remainder = chrLength % binSize
                start = int(np.floor(remainder/2)) 
                end = start + binSize
                # Calculate start and stop of all bins
                binEnd = np.arange(end, int(chrLength) + 1, binSize,
                    dtype = np.uint32)
                binStart = (binEnd - binSize)
            # Generate bin start and end sites for unequal width bins
            else:
                # Calculate numbe and size of bins
                binNo = np.ceil(float(chrLength) / binSize)
                excess = (binNo * binSize) - int(chrLength)
                largeBin = int(binSize - np.floor(excess / binNo))
                smallBin = int(binSize - np.ceil(excess / binNo))
                # Calculate bin numbers
                smallBinNo = int(excess % binNo)
                largeBinNo = int(binNo - smallBinNo)
                # Generate bins
                binWidth = np.array([largeBin] * largeBinNo + 
                    [smallBin] * smallBinNo)
                binEnd = np.cumsum(binWidth, dtype = np.uint32)
                binStart = (binEnd - binWidth)
            # Generate arrays containing bin index and chromosome name
            count += len(binStart)
            binCounts = np.zeros(len(binStart), dtype = np.uint32)
            # Store output data
            binDict[chrName] = {'start' : binStart, 'end' : binEnd, 'count' :
                binCounts}
        # Return dataframe
        return(binDict)
        
    def add_bin_count(self, binDict, chrom, position):
        ''' Return index of bin of specified chromosme position.
        
        Args:
            binDict (dict)- A dictionary of bin data generated using the
                self.create_bins function.
            chrom (str)- Chromsome name.
            position (int)- Position on chromosome. 0 index half-closed.
        
        Returns:
            binIndex (int)- Index of bin as defined within binDict.
        
        '''
        # Convert position format
        position = np.uint32(position)
        binLocation = binDict[chrom]['end'].searchsorted(
            position, side='right')
        # Extract and return bin index
        try:
            binStart = binDict[chrom]['start'][binLocation]
        except IndexError:
            binStart = position + 1
        # Extract bin name/number
        if binStart <= position:
            binDict[chrom]['count'][binLocation] += 1
    
    def count_bin_overlaps(
            self, binSize, binEqual, mapQ, overlap, rmDup=False, rmSec=False,
            rmSup=False
        ):
        ''' Counts overlaps between reads in a BAM file and genome bins.
        
        Args:
            binSize (int)- Maximum size of bins.
            binEqual (bool)- Whether bins should be equally sized.
            overlap (str)- Must be one of the following:
                'refstart' - 5' most portion of the read on the reference.
                'refend' - 3' most portion of the read on the reference.
                'readstart' - First sequenced base of the read.
                'readend' - Last sequence base of the read.
            rmDup (bool)- Remove duplicate reads
            rmSecond (bool)- Remove secondary reads.
        
        '''
        # Check arguments
        if not isinstance(rmDup, bool):
            raise TypeError('rmDup must be bool')
        if not isinstance(rmSec, bool):
            raise TypeError('rmSecond must be bool')
        if not isinstance(rmSup, bool):
            raise TypeError('rmSupplement must be bool')
        # Create filter flag
        flagFilter = 516
        if rmDup:
            flagFilter += 1024
        if rmSec:
            flagFilter += 256
        if rmSup:
            flagFilter += 2048
        # Create bin dictionary and output
        binDict = self.create_bins(binSize, binEqual)
        # Open bam and loop through
        inbam = pysam.AlignmentFile(self.bam)
        for read in inbam:
            # Skip unwanted and poorly mapped reads
            if read.flag & flagFilter:
                continue
            if read.mapping_quality < mapQ:
                continue
            # Extract chromosome name
            ref_name = read.reference_name
            ref_start = read.reference_start
            ref_end = read.reference_end - 1
            if overlap == 'refstart':
                self.add_bin_count(binDict, ref_name, ref_start)
            elif overlap == 'refend':
                self.add_bin_count(binDict, ref_name, ref_end)
            elif overlap == 'readstart':
                if read.is_reverse:
                    self.add_bin_count(binDict, ref_name, ref_end)
                else:
                    self.add_bin_count(binDict, ref_name, ref_start)
            elif overlap == 'readend':
                if read.is_reverse:
                    self.add_bin_count(binDict, ref_name, ref_start)
                else:
                    self.add_bin_count(binDict, ref_name, ref_end)
            else:
                raise ValueError('Unrecognised value for overlap')
        # Return data
        return(binDict)

class multiple_coverage(object):
    
    def __init__(self, bamList):
        ''' Function to initialise multiple_coverage object. Object
        functions to perform coverage calculations across BAM files.
        All BAM files must have reference sequences of the same name
        and length.
        
        Args:
            bamList (iter)- Iterable returning path to BAM file.
        
        Raises:
            ValueError - If reference sequences within BAM files
                do not have the same names and lengths.
        
        '''
        # Check bams contain identical reference sequence
        for count, bam in enumerate(bamList):
            bamCov = single_coverage(bam)
            if count:
                if bamCov.length != reference:
                    raise ValueError('Chromosomes in BAM files vary')
            else:
                reference = bamCov.length
        # Store bam list and lengths
        self.length = reference
        self.bamList = bamList
    
    def mean_coverage(
        self, intervals, map_quality=0, remove_dup=False
    ):
        ''' Function to return mean coverage across all intervals
        within matched BAM files.
        
        Args:
            intervals - Iterable of intervals where each element consists
                of a string and two integers specifying chromosome name and
                chromosome start and end, respectively.
            map_quality (int)- Minimum mapping quality for reads.
            remove_dup (bool)- Remove duplicate reads.
        
        Returns:
            outDF - Mean coverage of all intervals in all BAM files.
        
        '''
        # Check intervals
        for count, bam in enumerate(self.bamList):
            # Generate single coverage object
            bamCov = single_coverage(bam)
            # Check intervals for first bam and build output
            if count == 0:
                bamCov.check_intervals(intervals)
                names = ['{}:{}-{}'.format(*x) for x in intervals]
                outDF = pd.DataFrame(index=names, columns=self.bamList)
            # Add coverage to output
            outDF[bam] = bamCov.mean_coverage_each(
                intervals=intervals, map_quality=map_quality,
                remove_dup=remove_dup)
        # Return data
        return(outDF)
    
    def mean_coverage_cor(
        self, intervals, map_quality=0, remove_dup=False, method='pearson',
        min_coverage=0
    ):
        ''' Function to calculate correlation between mean coverage of
        intervals within BAM files.
        
        Args:
            intervals - Iterable of intervals where each element consists
                of a string and two integers specifying chromosome name and
                chromosome start and end, respectively.
            map_quality (int)- Minimum mapping quality for reads.
            remove_dup (bool)- Remove duplicate reads.
            method (str)- Method for correlation calculation.
            min_coverage (int)- Minimum coverage of interval, across all BAM
                files, for inclusion in correlation calculation.
        
        Returns:
            cor - Correlation matrix for bamFiles.
            meanCov - Mean coverage of all intervals in all BAM files.
        
        '''
        # Check arguments
        if not isinstance(return_filtered, bool):
            raise TypeError('return_filtered must be bool')
        methods = ('pearson', 'kendall', 'spearman')
        if method not in methods:
            raise ValueError('method but be one of: {}.'.format(
                ', '.join(methods)))
        # Extract mean coverage and remove low rows
        meanCov = self.mean_coverage(
            intervals=intervals, map_quality=map_quality,
            remove_dup=remove_dup)
        # Perform correlation between remaining points and return
        cor = filtCov.corr(method)
        return((cor, meanCov))
