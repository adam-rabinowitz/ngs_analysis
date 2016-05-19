'''bam2bedgraph.py

Usage:
    
    bam2bed.py <bam> <bed>

'''
# Import required modules
import pysam
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__, version='1.0')
# Open input bam file and output bed
bamFile = pysam.AlignmentFile(args['<bam>'])
bedFile = open(args['<bed>'], 'w')
# Create chromosome dictionary
chromDict = {}
for chrom in bamFile.references:
    chromDict[bamFile.gettid(chrom)] = chrom
# Create strand and read dictionary
strandDict = {False:'+', True:'-'}
readDict = {True:'/1', False:'/2'}
# Create output bed
for read in bamFile:
    if read.is_unmapped:
        continue
    bedFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(
        chromDict[read.reference_id],
        read.reference_start,
        read.reference_end,
        read.query_name + readDict[read.is_read1],
        read.mapping_quality,
        strandDict[read.is_reverse]
    ))
# Close input and output bed files
bamFile.close()
bedFile.close()
