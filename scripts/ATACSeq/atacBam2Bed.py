"""atacBam2Bed.py

Usage:
    
    atacBam2Bed.py <bamfile> <outprefix> [--minMapQ=<minMapQ>] [--size=<size>]
        [--sort]
    
    atacBam2Bed.py (-h | --help)
    
Options:
    
    --minMapQ=<minMapQ>  Minimum mapping quality of reads [default: 20].
    --size=<size>        Size of open region around insertion [default: 50].
    --sort               Sort bedfile. Requires bedtools to be in path.
    --help               Output this message.
    
"""
import pysam
import collections
import subprocess
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
    raise IOError('size argument must be divisible by 2')
# Process size argument
try:
    args['--minMapQ'] = int(args['--minMapQ'])
except ValueError:
    raise IOError('minMapQ argument must be an integer')
# Create file names and open input files
tempBedPath = args['<outprefix>'] + '.temp.bed'
finalBedPath = args['<outprefix>'] + '.bed'
tempBed = open(tempBedPath, 'w')
inBam = pysam.AlignmentFile(args['<bamfile>'])
# Create variables for processing BAM files
chrDict = pysamfunc.chr_dict(inBam)
strDict = {True: '-', False: '+'}
counter = collections.defaultdict(int)
# Loop through BAM file and create BED file
for read in inBam:
    # Count total and extract flag
    counter['total'] += 1
    flag = read.flag
    # Count and skip unmapped reads
    if flag & 4:
        counter['unmapped'] += 1
        continue
    # Count and skip poorly mapped reads
    if read.mapping_quality < args['--minMapQ']:
        counter['poormap'] += 1
        continue
    # Count and skip reads with unmapped mate
    if flag & 8:
        counter['nomate'] += 1
        continue
    # Count and skip improper pairs
    if not flag & 2:
        counter['improper']
    # Count and skip duplicates
    if flag & 1024:
        counter['duplicate'] += 1
        continue
    # Count and skip secondary alignments
    if flag & 256:
        counter['secondary'] += 1
        continue
    # Count and skip secondary alignments
    if flag & 2048:
        counter['supplement'] += 1
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
    # Count and write accepted reads
    counter['accepted'] += 1
    positionTuple = pysamfunc.startInterval(read, args['--size'], chrDict)
    tempBed.write('\t'.join(map(str, positionTuple)) + '\n')
# Close input file
tempBed.close()
inBam.close()
# Open log file and write
with open(args['<outprefix>'] + '.log', 'w') as logfile:
    logfile.write('Total: {}\n'.format(counter['total']))
    logfile.write('Unmapped: {}\n'.format(counter['unmapped']))
    logfile.write('Map quality <{}: {}\n'.format(
        args['--minMapQ'], counter['poormap']))
    logfile.write('Mate Unmapped: {}\n'.format(counter['nomate']))
    logfile.write('Improper pair: {}\n'.format(counter['improper']))
    logfile.write('Duplicate: {}\n'.format(counter['duplicate']))
    logfile.write('Secondary: {}\n'.format(counter['secondary']))
    logfile.write('Supplementary: {}\n'.format(counter['supplement']))
    logfile.write('Start clipped: {}\n'.format(counter['clipped']))
    logfile.write('Accepted: {}\n'.format(counter['accepted']))
# Create and submit sort command
if args['--sort']:
    sortCommand = 'bedtools sort -i {} > {} && rm {}'.format(
        tempBedPath, finalBedPath, tempBedPath)
else:
    sortCommand = 'mv {} {}'.format(tempBedPath, finalBedPath)
subprocess.check_call(sortCommand, shell=True)
