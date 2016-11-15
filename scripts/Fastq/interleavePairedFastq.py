'''interleavePairedFastq.py

Usage:
    
    interlevePairedFastq.py <inFastq1> <inFastq2> <outFastq>
    
'''
# Import arguments
import os
from ngs_python.fastq import fastqIO
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Create input fastq objects objects
pf = fastqIO.parseFastq(
    fastq1=args['<inFastq1>'], fastq2=args['<inFastq2>'])
# Interleave fastq files
count = pf.interleave_reads(args['<outFastq>'])
print(count)
