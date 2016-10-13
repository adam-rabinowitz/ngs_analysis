'''sampleFastq.py

Usage:
    
    sampleFastq.py single <inFastq1> <outFastq1> <number>
    sampleFastq.py pair <inFastq1> <inFastq2> <outFastq1> <outFastq2>
        <number>
    
'''
# Import arguments
import os
from ngs_python.fastq import fastqIO
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<number>'] = int(args['<number>'])
# Create fastq object and count and check reads
if args['single']:
    fastq = fastqIO.readFastq(args['<inFastq1>'])
elif args['pair']:
    fastq = fastqIO.readFastq(args['<inFastq1>'], args['<inFastq2>'])
count = fastq.check_names()
print(count)
if count < args['<number>']:
    raise ValueError('Number exceeds reads in FASTQ file')
# Sample FASTQ files
if args['pair']:
    fastq.sample_reads(
        number=args['<number>'], sample=count, outFastq1=args['<outFastq1>'],
        outFastq2 = ['<outFastq2>'])
elif args['single']:
    fastq.sample_reads(
        number=args['<number>'], sample=count, outFastq1=args['<outFastq1>'])
