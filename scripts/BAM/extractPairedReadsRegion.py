'''extractPairedReadsRegion.py
    
Usage:
    extractPairedReadsRegion.py <region> <inbam> <outbam>
        [--threads=<th>] [--memory=<mm>] [--nosort]

Options:
    --threads=<th>  Number of threads to use in sort [default: 1].
    --memory=<mm>   GB of memory per sort thread [default: 6].
    --nosort        Input BAM is already sorted by name.

'''
# Load required modules
from general_python import docopt
from ngs_python.bam.samtools import sort as bamsort
import os
import pysam
import subprocess
# Process arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--threads'] = int(args['--threads'])
if args['--threads'] < 1:
    raise ValueError('Must be at least 1 thread')
args['--memory'] = int(args['--memory'])
if args['--memory'] < 1:
    raise ValueError('Must be at least 1 GB of memory')
args['<chrom>'], interval = args['<region>'].split(':')
args['<start>'], args['<end>'] = interval.split('-')
args['<start>'] = int(args['<start>']) - 1
args['<end>'] = int(args['<end>'])
if not args['<inbam>'].endswith('.bam'):
    raise IOError('Unexpected input file name')
args['<inbam>'] = os.path.abspath(args['<inbam>'])
if not args['<outbam>'].endswith('.bam'):
    raise IOError('Unexpected output file name')
args['<outbam>'] = os.path.abspath(args['<outbam>'])
# Name sort input BAM file if required
if args['--nosort']:
    args['<nsortbam>'] = args['<inbam>']
else:
    args['<nsortbam>'] = args['<outbam>'][:-4] + '.nsort.bam'
    inputSortCommand = bamsort(
        inFile=args['<inbam>'], outFile=args['<nsortbam>'], name=True,
        threads=args['--threads'], memory=args['--memory'], delete=False
    )
    subprocess.check_call(inputSortCommand, shell=True)
# Create temporary output bam file and extract tid
args['<filtbam>'] = args['<outbam>'][:-4] + '.filtered.bam'
sortbam = pysam.AlignmentFile(args['<nsortbam>'], 'rb')
args['<tid>'] = sortbam.get_tid(args['<chrom>'])
filtbam = pysam.AlignmentFile(args['<filtbam>'], 'wb', sortbam)
# Extract reads from input BAM
currentName = ''
readList = []
while True:
    # Extract read data
    try:
        nextRead = sortbam.next()
        nextName = nextRead.query_name
    except StopIteration:
        nextName = None
    # Process reads of the same name
    if nextName != currentName:
        # Determine if one read maps to region of interest
        passed = False
        for read in readList:
            if read.is_unmapped:
                continue
            elif read.reference_id != args['<tid>']:
                continue
            elif read.reference_start >= args['<end>']:
                continue
            elif read.reference_end <= args['<start>']:
                continue
            else:
                passed = True
                break
        # Write mapped reads to file
        if passed:
            for read in readList:
                filtbam.write(read)
        # Reset read list and current name
        if nextName == None:
            break
        else:
            currentName = nextName
            readList = [nextRead]
    # Create list of reads of the same name
    else:
        readList.append(nextRead)
# Close BAM files
sortbam.close()
filtbam.close()
# Delete temporary name sorted bam
if not args['--nosort']:
    os.remove(args['<nsortbam>'])
# Sort output bam file
outputSortCommand = bamsort(
    inFile=args['<filtbam>'], outFile=args['<outbam>'], name=False,
    threads=args['--threads'], memory=args['--memory'], delete=True
)
subprocess.check_call(outputSortCommand, shell=True)
