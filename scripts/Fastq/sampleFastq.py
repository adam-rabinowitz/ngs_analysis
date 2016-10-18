'''sampleFastq.py

Usage:
    
    sampleFastq.py single <number> <inFastq1> <outFastq1>
    sampleFastq.py pair <number> <inFastq1> <inFastq2> <outFastq1> <outFastq2>
    
'''
# Import arguments
import os
from ngs_python.fastq import fastqIO
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<number>'] = int(args['<number>'])
if args['<number>'] < 1:
    raise ValueError('Number must be positive integer')
# Create input fastq objects objects
parseIn = fastqIO.parseFastq(
    fastq1=args['<inFastq1>'], fastq2=args['<inFastq2>'])
countIn = parseIn.check_names()
print '{} initial reads'.format(countIn)
if countIn < args['<number>']:
    raise ValueError('Number exceeds reads in FASTQ file')
# Sample FASTQ files
parseIn.sample_reads(
    number=args['<number>'], sample=countIn, outFastq1=args['<outFastq1>'],
    outFastq2 = args['<outFastq2>'])
# Check output reads
parseOut = fastqIO.parseFastq(args['<outFastq1>'], args['<outFastq2>'])
countOut = parseOut.check_names()
print('{} final reads'.format(countOut))
if countOut != args['<number>']:
    raise ValueError('Number does not equal output count')
