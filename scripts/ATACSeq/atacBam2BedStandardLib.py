"""atacBam2Bed.py

Usage:
    
    atacBam2Bed.py <bamfile> <size> <outbed> [--minMapQ=<minMapQ>]
        [--rmDup]
    
    atacBam2Bed.py (-h | --help)
    
Options:
    
    --minMapQ=<minMapQ>  Minimum mapping quality of reads [default: 20]
    --size=<size>        Size of open region around insertion [default: 50]
    --rmDup              Flag to remove duplicate reads
    --help               Output this message
    
"""
import argparse
import collections
import pysam
# Process arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    'bamfile', type=str, nargs='+', help='input bam file'
)
parser.add_argument(
    'outbed', type=str, nargs='+', help='output bed file'
)
parser.add_argument(
    'size', type=int, nargs='+', help='size of window'
)
parser.add_argument(
    'minMapQ', type=int, nargs='+', help='minimum mapping quality'
)
parser.add_argument(
    '--rmDup', type=bool, action='store', help='Remove duplicate reads'
)
                                        const=sum, default=max,
                                                            help='sum the integers (default: find the max)')

# Function to extract start interval data from read
def startInterval(read, intervalSize, chrDict):
    ''' Function returns interval around the start of a read. If the
    read is aligned to forward strand then the start is defined as the
    end of the read aligned to the most 5' portion of the reference.
    If the read is aligned to the reverse strand then the start is defined
    as end of the read aligned to the most 3' portion of the reference.
    The function takes the following 3 arguments:
    
    1)  read - Pysam AlignedSegment
    2)  intervalSize - Size of the returned interval
    3)  chrDict - A dictionary where the keys are the target IDs (tid)
        of the chromosomes and the values is a tuple of the chromsome
        name and length.
    
    Function returns a tuple of the region in BED format; using zero
    based indexing and open ended interval and where name is read name
    and score is read mapping quality.
    
    '''
    # Extract flag
    flag = read.flag
    # Find strand of read and extract read start position
    if flag & 16:
        strand = '-'
        site = read.reference_end
    else:
        strand = '+'
        site = read.reference_start
    # Calculate start and end of interval
    halfSize = int(intervalSize / 2)
    start = site - halfSize
    end = site + halfSize
    # Extract chromosome data and adjust start and end
    chrom, length = chrDict[read.reference_id]
    start = max(0, start)
    end = min(length, end)
    # Extract mapping quality
    mapq = read.mapping_quality
    # Extract read name
    if flag & 64:
        name = read.query_name + '/1'
    elif flag & 128:
        name = read.query_name + '/2'
    else:
        name = read.query_name
    # Return data
    return((chrom, start, end, name, mapq, strand))

# 
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
chrDict = dict([(inBam.get_tid(x),(x,y)) for x,y in zip(
    inBam.references, inBam.lengths)])
# Loop through BAM file and create BED file
for read in inBam:
    # Count total and extract flag
    counter['total'] += 1
    flag = read.flag
    # Count and skip duplicates
    if args['--rmDup'] and flag & 1024:
        counter['duplicate'] += 1
        continue
    # Count and skip secondary alignments
    if flag & 256:
        counter['secondary'] += 1
        continue
    # Count and skip unmapped reads
    if flag & 4:
        counter['unmapped'] += 1
        continue
    # Count and skip reads with unmapped mate
    if flag & 8:
        counter['nomate'] += 1
        continue
    # Count and skip poorly mapped reads
    if read.mapping_quality < args['--minMapQ']:
        counter['poormap'] += 1
        continue
    # Count and skip reads with clipping at read start
    if flag & 16:
        if read.cigartuples[-1][0] in [4, 5]:
            counter['clipped'] += 1
            continue
    else:
        if read.cigartuples[0][0] in [4, 5]:
            counter['clipped'] += 1
            continue
    # Count accepted reads
    counter['accepted'] += 1
    # Extract bed data from read
    regionTuple = startInterval(read, args['--size'], chrDict)
    # Save regions to file
    regionTuple = '\t'.join(map(str, regionTuple))
    outBed.write('%s\n' %(regionTuple))
# Close files
outBed.close()
inBam.close()
# Print metric
print 'Total reads: %s' %(counter['total'])
print 'Duplicates: %s' %(counter['duplicate'])
print 'Secondary: %s' %(counter['secondary'])
print 'Unmapped: %s' %(counter['unmapped'])
print 'Mate unmapped: %s' %(counter['nomate'])
print 'Poorly mapped: %s' %(counter['poormap'])
print 'Start clipped: %s' %(counter['clipped'])
print 'Accepted: %s' %(counter['accepted'])
