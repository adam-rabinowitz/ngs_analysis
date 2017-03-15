import collections
import os
import pysam

class PysamBAM(object):
    
    def __init__(self, bam, fasta = None):
        ''' Initiates BAMAnalysis object. Extracts names and lengths of
        reference sequences for supplied BAM and FASTA files and check for
        consistency.
        
        Args:
            bam - Path to BAM file or an iterable of BAM file paths.
            fasta - Path to reference FASTA file for BAM file(s).
        
        '''
        # Extract absolute paths for bam files and check
        bampaths = []
        if isinstance(bam, str):
            bampaths.append(bam)
        else:
            for b in bam:
                bampaths.append(b)
        bampaths = [os.path.abspath(bp) for bp in bampaths]
        for bp in bampaths:
            if not os.path.isfile(bp):
                raise IOError('Could not find file {}'.format(bp))
        self.bampaths = bampaths
        # Check all bam files for index and consistency
        for count, bp in enumerate(bampaths):
            with pysam.AlignmentFile(bp, 'rb') as inbam:
                # Check index
                if not inbam.check_index():
                    raise IOError('No index for file {}'.format(bp))
                # Check reference sequence consistency
                lengths = collections.OrderedDict(
                    zip(inbam.references, inbam.lengths))
                if count:
                    if lengths != self.lengths:
                        raise ValueError('BAM file sequences differ')
                else:
                    self.lengths = lengths
        # Extract and check names and lengths of FASTA sequences
        if fasta is not None:
            fastapath = os.path.abspath(fasta)
            with pysam.FastaFile(fastapath) as infasta:
                lengths = collections.OrderedDict(
                    zip(infasta.references, infasta.lengths))
                if lengths != self.lengths:
                    raise ValueError('FASTA file sequences differ')
            self.fasta = fastapath
        else:
            self.fasta = None
    
    def __enter__(self):
        return(self)
    
    def __exit__(self, exc_type, exc_value, traceback):
        pass
    
    def _flagfilter(
            self, unmapped = False, qcfail = False, duplicate = False,
            secondary = False, supplementary = False, properpair = True
        ):
        ''' Creates flag to filter bam reads.
        
        Args:
            unmapped (bool)- Return unmapped reads.
            qcfail (bool)- Return qc failed reads.
            duplicate (bool)- Return duplicate alignments.
            secondary (bool)- Return secondary alignments.
            supplementary (bool)- Return secondary alignments.
        
        Returns:
            filter_flag (int)- Bitwise flag for filtering BAM files.
        
        '''
        # Check arguments
        if not isinstance(unmapped, bool):
            raise TypeError('unmapped must be bool')
        if not isinstance(qcfail, bool):
            raise TypeError('qcfail must be bool')
        if not isinstance(duplicate, bool):
            raise TypeError('duplicate must be bool')
        if not isinstance(secondary, bool):
            raise TypeError('secondary must be bool')
        if not isinstance(supplementary, bool):
            raise TypeError('supplementary must be bool')
        if not isinstance(properpair, bool):
            raise TypeError('properpair must be bool')
        # Build negative filter
        neg_filter = 0
        if not unmapped:
            neg_filter += 4
        if not secondary:
            neg_filter += 256
        if not qcfail:
            neg_filter += 512
        if not duplicate:
            neg_filter += 1024
        if not supplementary:
            neg_filter += 2048
        # Build positive filter
        pos_filter = 0
        if properpair:
            pos_filter += 2
        # Create and return filter function
        def filterfunc(flag):
            if flag & neg_filter:
                return(True)
            elif (flag & pos_filter) != pos_filter:
                return(True)
            else:
                return(False)
        return(filterfunc)
        
    def _getreads(
            self, bam, filterfunc, chrom = None, start = None, end = None,
            flag = 3844, mapq = None
        ):
        ''' Creates generator returning aligned reads from BAM file.
        
        Args:
            bam (str)- A pysam.AlignmentFile object.
            chrom (str)- Chromosome from which to extract reads.
            start (int)- Start of region from which to extract reads.
            end (int)- End of region from which to extract reads.
            flag (int)- Bitwise flag for filtering reads.
            mapq (int)- Minimum read map quality.
        
        '''
        # Check arguments
        if mapq is not None:
            if not isinstance(mapq, int):
                raise TypeError('mapq must be integer')
            if mapq < 0:
                raise ValueError('mapq must be non-negative')
            if not flag & 4:
                raise ValueError('Must filter unmapped reads when mapq set')
        # Create generator loop
        for read in bam.fetch(chrom, start - 1, end):
            if read.flag & flag:
                continue
            if mapq and read.mapping_quality < mapq:
                continue
            yield(read)
    
    def _getpileup(
            self, bam, filterfunc, chrom, pos, mapq = None
        ):
        ''' Creates generator returning aligned reads from BAM file.
        
        Args:
            bam (str)- A pysam.AlignmentFile object.
            filterfunc (function)- Function to filter read flags. Generated
                by self._flagfilter.
            chrom (str)- Chromosome from which to extract reads.
            start (int)- Start of region from which to extract reads.
            end (int)- End of region from which to extract reads.
            flag (int)- Bitwise flag for filtering reads.
            mapq (int)- Minimum read map quality.
        
        '''
        # Check arguments
        if mapq is not None:
            if not isinstance(mapq, int):
                raise TypeError('mapq must be integer')
            if mapq < 0:
                raise ValueError('mapq must be non-negative')
            if not flag & 4:
                raise ValueError('Must filter unmapped reads when mapq set')
        # Adjust position
        pos -= 1
        # Create pileup generator
        for column in bam.pileup(chrom, pos, pos + 1, stepper='all'):
            # Loop through reads in desired column
            if column.pos == pos:
                for read in column.pileups:
                    # Skip reads without base calls
                    if read.is_del or read.is_refskip:
                        continue
                    # Extract query data
                    alignment = read.alignment
                    query_position = read.query_position
                    # Skip unwanted reads
                    if filterfunc(alignment.flag):
                        continue
                    # Check mapping quality
                    if mapq and alignment.mapping_quality < mapq:
                        continue
                    # Return read data
                    yield(alignment, query_position)
    
    def check_intervals(self, intervals):
        ''' Checks supplied intervals for validity.
        
        Args:
            intervals - Iterator returning chromosome name, start position
                and end position of intervals.
        
        Raises:
            TypeError - If name, start and end are not a string, integer and
                integer, respectively.
            IndexError - If name is not a valid reference sequence.
            ValueError - If start and end values are invalid relative to each
                other or reference length.
        
        '''
        # Iterate through intervals and check values and types
        for chrom, start, end in intervals:
            # Check variable types
            if not isinstance(chrom, str):
                raise TypeError('Interval chromosomess must be strings')
            if not isinstance(start, int):
                raise TypeError('Interval starts must be integers')
            if not isinstance(end, int):
                raise TypeError('Interval ends must be integers')
            # Check variable values
            chromLength = self.lengths[chrom]
            if start > end:
                raise ValueError('Interval start > interval end')
            if start < 1:
                raise ValueError('Interval start not on chromosome')
            if end > chromLength:
                raise ValueError('Interval end not on chromosome')

