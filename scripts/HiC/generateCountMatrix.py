"""generateCountMatrix.py

Usage:
    
    generateCountMatrix.py bed <bedfile> <infile> <outfile>
        [--threads=<threads>]
    
    generateCountMatrix.py nobed <chrfile> <binsize> <infile>
        <outfile> [--threads=<threads>] [--equal]
    
    generateCountMatrix.py (-h | --help)
    
Options:
    
    --threads=<threads>  Number of threads [default: 1]
    --equal              Bins should be equally sized
    --help               Output this message
    
"""
# Import required modules
import os
import re
import numpy as np
from ngs_python.structure import interactionMatrix
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Check numerical arguments
args['--threads'] = int(args['--threads'])
if args['nobed']:
    args['<binsize>'] = int(args['<binsize>'])
# Check input files
if not os.path.isfile(args['<infile>']):
    raise IOError('{} not found'.format(args['<infile>']))
if args['bed']:
    if not os.path.isfile(args['<bedfile>']):
        raise IOError('{} not found'.format(args['<bedfile>']))
else:
    if not os.path.isfile(args['<chrfile>']):
        raise IOError('{} not found'.format(args['<chrfile>']))
# Extract and print parameters to create bins
if args['bed']:
    binData = args['<bedfile>']
    print '\nParameters:\n  %s\n' %(
        'bed file provided',
    )
else:
    binData = (args['<chrfile>'], args['<binsize>'], args['--equal'])
    print '\nParameters:\n  %s\n  %s\n' %(
        'max bin size: %s' %(args['<binsize>']),
        'bin size equal: %s' %(args['--equal'])
    )
# Create bin object and save bed
genomeBins = interactionMatrix.genomeBin(binData)
# Create interaction matrix and save to file
countMatrix, logArray = genomeBins.generateMatrix(args['<infile>'],
    args['--threads'])
np.savetxt(args['<outfile>'], countMatrix, fmt = '%s', delimiter = '\t',
    header = '\t'.join(genomeBins.binNames), comments = '')
# Print interaction data
print 'Interaction Data:\n  %s\n  %s\n  %s\n  %s\n' %(
    'total: %s' %(logArray[0]),
    'accepted: %s' %(logArray[3]),
    'no chromosome: %s' %(logArray[1]),
    'no bin: %s' %(logArray[2])
)
