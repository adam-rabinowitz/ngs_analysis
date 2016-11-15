'''extractChrSizes.py

Usage:
    
    extractChrSizes.py <fasta> <out>
    
'''
from Bio import SeqIO
from general_python import docopt
args = docopt.docopt(__doc__,version = 'v1')
with open(args['<out>'], 'w') as outFile:
    for sequence in SeqIO.parse(open(args['<fasta>']), 'fasta'):
        outFile.write('{}\t{}\n'.format(sequence.id, len(sequence)))
