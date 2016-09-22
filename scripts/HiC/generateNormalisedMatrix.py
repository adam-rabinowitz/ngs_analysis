"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <regions> <mincount> <outdir> <infiles>...
        [--path=<paths>] [--threads=<threads>] [--memory=<memory>]
        [--iter=<iter>]
    
    generateNormalisedMatrix.py (-h | --help)
    
Options:
    
    --path=<path>        Path to ic_mes executable [default: ic_mes]
    --threads=<threads>  Number of threads [default: 1]
    --memory=<memory>    Memory (MB) per thread [default: 2000]
    --iter=<iter>        Iterations of ic_mes [default: 1000]
    --help               Output this message
    
"""
# Import required modules
from ngs_python.structure import interactionMatrix
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__, version = '1.0')
# Check numeric arguments
args['<mincount>'] = int(args['<mincount>'])
args['--threads'] = int(args['--threads'])
args['--memory'] = int(args['--memory'])
args['--iter'] = int(args['--iter'])
# Perform processing
normMatrices = interactionMatrix.normaliseCountMatrices(
    args['<infiles>'], args['<regions>'])
# Find low bins
normMatrices.extractLowBins(args['<mincount>'], args['--threads']) 
# Create sub-matrices
normMatrices.saveSubMatrices(args['<outdir>'], args['--threads'])
# Normalise sub matrices
normMatrices.normaliseSubMatrices(args['--path'], args['--threads'])
