"""extractQuantileDistance.py

Usage:
    
    extractQuantileDistance.py <quantile> <outfile> <infiles>...
        [--threads=<threads>]
    
    extractQuantileDistance.py (-h | --help)
    
Options:
    
    --threads=<threads>  Number of threads [default: 4]
    
"""
# Import modules
import os
from general_python import docopt
from ngs_python.structure import analyseInteraction
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<quantile>'] = float(args['<quantile>'])
if not 1 > args['<quantile>'] > 0:
    raise ValueError('quantile must be between 0 and 1')
args['--threads'] = int(args['--threads'])
if not 13 > args['--threads'] > 0:
    raise ValueError('threads must be >= 0 and <= 12')
outDir = os.path.dirname(args['<outfile>'])
if not os.path.isdir(outDir):
    raise ValueError('Could not find path for outfile')
# Create 
ai = analyseInteraction.analyse_interaction(args['<infiles>'])
quantileData = ai.calc_quantile(args['<quantile>'], args['--threads'])
quantileData.to_csv(args['<outfile>'], index=False, sep='\t')
