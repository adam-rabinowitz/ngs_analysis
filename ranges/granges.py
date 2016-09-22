import pysam
import pandas as pd
import collections
import copy

class granges:

    def __init__(self, lengths = None):
        # Read in chromosome lengths
        if lengths.endswith('.fa') or lengths.endswith('.fasta'):
            in_file = pysam.FastaFile(lengths)
            self.lengths = collections.OrderedDict(
                zip(in_file.references, in_file.lengths)
            )
            in_file.close()
        elif lengths.endswith('.bam'):
            in_file = pysam.BamFile(lengths)
            self.lengths = collections.OrderedDict(
                zip(in_file.references, in_file.lengths)
            )
            in_file.close()
        elif lengths is not None:
            self.lengths = collections.OrderedDict()
            with open(lengths) as in_file:
                for line in in_file:
                    chrom, length = line.strip().split('\s')[0:2]
                    self.lengths[chrom] = int(length)
        else:
            self.lengths = None
        if self.lengths:
            self.intervals = collections.OrderedDict()
            for chrom in self.lengths:
                self.intervals[chrom] = pd.DataFrame()
    
    def add_bed_intervals(self, bed):
        # Read in bed data and add column names
        bed_data = pd.DataFrame.from_csv(path = bed, sep = '\t',
            header = None, index_col = False
        )
        if bed_data.shape[1] == 3:
            bed_data.columns = ['chr', 'start', 'end']
            bed_data['strand'] = '*'
        elif bed_data.shape[1] == 6:
            bed_data.columns = ['chr', 'start', 'end', 'name', 'score',
                'strand']
        else:
            raise ValueError('bed files must have 3 or 6 columns')
        if (bed_data['start'] > bed_data['end']).any():
            raise ValueError('start greater than end')
        # Add bed data to interval list and sort
        bed_data = bed_data.groupby('chr')
        for chrom, chrom_data in bed_data:
            if self.lengths:
                if chrom not in self.lengths:
                    raise ValueError('chromosome not recognised')
                if (chrom_data['end'].max() > self.lengths[chrom]
                    or chrom_data['start'].min() < 0):
                    raise ValueError('interval not within chromosome')
            chrom_data.drop('chr', axis=1, inplace = True)
            self.intervals[chrom] = pd.concat(
                (self.intervals[chrom], chrom_data),
                axis = 0
            )
        self.sort_intervals()
    
    def sort_intervals(self):
        count = 0
        for chrom in self.intervals:
            nrow = self.intervals[chrom].shape[0]
            if nrow == 0:
                continue
            self.intervals[chrom].sort_values(['start', 'end'], inplace=True)
            self.intervals[chrom].index = range(count, count + nrow)
            count += nrow
    
    def merge_overlaps(self, intervals):
        # Group intervals into overlapping groups
        intervals = intervals[['start','end']]
        intervals['seperate'] = intervals['start'] >= intervals['end'].shift(1)
        intervals['group'] = intervals['seperate'].cumsum() - 1
        intervals = intervals.groupby('group')
        # Creeate and return output data frame
        output = pd.DataFrame()
        output['start'] = intervals['start'].min()
        output['end'] = intervals['end'].max()
        return(output)
    
    def reduce_intervals(self, merge = False, ignore_strand = False):
        for chrom in self.intervals:
            nrow = self.intervals[chrom].shape[0]
            if nrow == 0:
                continue
            # Merge all overlapping intervals
            if ignore_strand:
            
            # Merge intervals on one strand
            else:
                x['overlap'] = x['start'] < x['end'].shift(1)
                x['group'] = (~x['overlap']).cumsum()
                x = x.groupby('group')
                out = 



       
    def extract_overlaps(self, other):
        output = copy.deepcopy(self)
        for chrom in self.intervals:
            # Extract chromosome data for both object or skip
            try:
                self_data = self.intervals[chrom]
                other_data = other.intervals[chrom]
            except KeyError:
                continue
            try:
                self_index = self_data.index[0]
                other_index = other_data.index[0]
            except IndexError:
                continue
            # Find overlaps or advance indices
            output_indices = []
            output_start = []
            output_end = []
            while True:
                try:
                    self_start, self_end = self_data.loc[
                        self_index, ['start', 'end']
                    ]
                    other_start, other_end = other_data.loc[
                        other_index, ['start', 'end']
                    ]
                except KeyError:
                    break
                if self_end <= other_start:
                    self_index += 1
                elif other_end <= self_start:
                    other_index += 1
                elif self_end > other_start:
                    output_indices.append(self_index)
                    output_start.append(max(self_start, other_start))
                    output_end.append(min(self_end, other_end))
                    self_index += 1
            if len(output_indices) > 0:
                output.intervals[chrom] = output.intervals[chrom].loc[
                    output_indices
                ]
                output.intervals[chrom]['start'] = output_start
                output.intervals[chrom]['end'] = output_end
        return(output)

test1 = granges('/farm/scratch/rs-bio-lif/rabino01/Elza/genome/mm10.fa')
test1.add_bed_intervals('/farm/home/rabino01/test1.bed')
test2 = granges('/farm/scratch/rs-bio-lif/rabino01/Elza/genome/mm10.fa')
test2.add_bed_intervals('/farm/home/rabino01/test2.bed')
overlap = test1.extract_overlaps(test2)
print(overlap.intervals['chr1'])
print(test1.intervals['chr1'])
