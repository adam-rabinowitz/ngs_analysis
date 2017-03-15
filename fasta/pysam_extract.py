import pysam

class fasta(object):
    
    def __init__(self, fasta):
        self.fasta = pysam.FastaFile(fasta)
        self.lengths = dict([(x,y) for x,y in zip(
            self.fasta.references, self.fasta.lengths)])
    
    def extract(self, chrom, start, end, upper=False):
        ''' Function to extract sequence of specified region.
        
        Args:
            chrom (str)- Name of chromosome.
            start (int)- Start of regions.
            end (int)- End of sequence.
        
        Returns:
            seq (str)- String of sequence.
        
        '''
        # Check arguments
        if end > self.lengths[chrom]:
            raise ValueError('End is greater than chromosome length')
        # Extract and return sequence
        seq = self.fasta.fetch(chrom, start, end)
        return(seq)
    
    def extract_fasta(self, outfasta, regions, upper=False):
        ''' Function to create FASTA file of extracted regions
        
        Args:
            outfasta (str)- Name of output FASTA file.
            regions - Iterable returning chrom, start and end values.
            upper (bool)- Convert sequences to upper case.
        
        '''
        # Open and create output file
        with open(outfasta, 'w') as outfile:
            for chrom, start, end in regions:
                name = '>{}:{}-{}\n'.format(chrom, start + 1, end)
                outfile.write(name)
                seq = self.extract(chrom, start, end)
                outfile.write(seq)
                outfile.write('\n')
    
    def extract_fastq(self, outfastq, regions, quality=30):
        ''' Function to create FASTQ file of extracted regions
        
        Args:
            outfasta (str)- Name of output FASTA file.
            regions - Iterable returning chrom, start and end values.
            upper (bool)- Convert sequences to upper case.
        
        '''
        # Open and create output file
        with open(outfastq, 'w') as outfile:
            for chrom, start, end in regions:
                name = '@{}:{}-{}'.format(chrom, start + 1, end)
                seq = self.extract(chrom, start, end)
                quality = chr(quality + 33) * (end - start)
                fastqEntry = '{}\n{}\n+\n{}\n'.format(
                    name, seq, quality)
                outfile.write(qualityString + '\n')

    def tile_region(
            self, chrom, start, end, size, overlap=0, complete=True
        ):
        ''' Function to create list of inervals tiling a region
        
        Args:
            chrom (str)- Name of chromosome.
            start (int)- Start of interval to tile.
            end (int)- End of interval to tile.
            size (int)- Size of each interval.
            overlap (int)- Overlap between intervals.
        
        Returns:
            regions (list)- A list of tuples of the region.
        
        '''
        # Check arguments
        if not all([isinstance(x, int) for x in (start, end, size, overlap)]):
            raise TypeError('start, end, size and overlap must be integers')
        if not all([x > -1 for x in start, end, size]):
            raise ValueError('start, end and size must be non-negative')
        if start >= end:
            raise ValueError('Start must be less than end')
        if end > self.lengths[chrom]:
            raise ValueError('End is greater than chromosome length')
        # Create intervals
        startList = range(start, end - overlap, size - overlap)
        if (startList[-1] + size) != end and complete:
            raise ValueError('Supplied parameters do not span region')
        # Create and return output
        outList = [(chrom, x, x + size) for x in startList]
        return(outList)
