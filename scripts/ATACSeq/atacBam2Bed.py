"""atacBam2Bed.py

Usage:
    
    atacBam2Bed.py <bamfile> <outbed> [--minMapQ=<minMapQ>] [--size=<size>]
        [--rmDup=<rmDup>]
    
    atacBam2Bed.py (-h | --help)
    
Options:
    
    --minMapQ=<minMapQ>  Minimum mapping quality of reads [default: 20]
    --size=<size>        Size of open region around insertion [default: 50]
    --rmDup=<rmDup>      Flag to remove duplicate reads
    --help               Output this message
    
"""
import pysam
import collections
from general_python import docopt
from ngs_python.bam import pysamfunc
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Process size argument
try:
    args['--size'] = int(args['--size'])
except ValueError:
    raise IOError('size argument must be an integer divisible by 2')
if args['--size'] % 2:
    raise IOError('size argument must be an integer divisible by 2')
# Process size argument
try:
    args['--minMapQ'] = int(args['--minMapQ'])
except ValueError:
    raise IOError('minMapQ argument must be an integer')
# Create counter to store processing metrics
counter = collections.defaultdict(int)
# Open input and output files
outBed = open(args['<outbed>'], 'w')
inBam = pysam.AlignmentFile(args['<bamfile>'])
# Create dictionary contianing chromosome names and lengths
chrDict = pysamfunc.createChrDict(inBam)
# Loop through BAM file and create BED file
for read in inBam:
    # Count total
    counter['total'] += 1
    # Count and skip duplicates
    if args['--rmDup'] and read.flag & 1024:
        counter['duplicate'] += 1
        continue
    # Count and skip poorly mapped reads
    if read.mapping_quality < args['--minMapQ']:
        counter['poormap'] += 1
        continue
    # Count and skip reads with unmapped mate
    if read.flag & 8:
        counter['nomate'] += 1
        continue
    # Count and skip reads with clipping at read start
    if read.flag & 16:
        if read.cigartuples[-1][0] in [4, 5]:
            counter['clipped'] += 1
            continue
    else:
        if read.cigartuples[0][0] in [4, 5]:
            counter['clipped'] += 1
            continue
    # Count and create accepted reads
    counter['accepted'] += 1
    regionTuple = pysamfunc.startInterval(read, args['--size'], chrDict)
    regionTuple = '\t'.join(map(str, regionTuple))
    outBed.write('%s\n' %(regionTuple))
# Close files
outBed.close()
inBam.close()
# Print metric
print 'Total reads: %s' %(counter['total'])
print 'Duplicates: %s' %(counter['duplicate'])
print 'Poorly mapped: %s' %(counter['poormap'])
print 'Mate unmapped: %s' %(counter['nomate'])
print 'Start clipped: %s' %(counter['clipped'])
print 'Accepted: %s' %(counter['accepted'])
