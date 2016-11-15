"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <regions> <mincount> <outdir> <infiles>...
        [--threads=<threads>] [--iter=<iter>] [--keepdiag]
    
    generateNormalisedMatrix.py (-h | --help)
    
Options:
    
    --threads=<threads>  Number of threads [default: 1]
    --iter=<iter>        Iterations of ic_mes [default: 1000]
    --keepdiag           Keep diagonal of matrix, otherwise set to zero.
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
args['--iter'] = int(args['--iter'])
# Perform processing
normMatrices = interactionMatrix.normaliseCountMatrices(
    matrixList=args['<infiles>'], regionFile=args['<regions>'])
# Normalise sub matrices
normMatrices.iceNormalisation(
    outDir=args['<outdir>'], minCount=args['<mincount>'],
    threads=args['--threads'], max_iter=args['--iter'],
    rmDiag=not(args['--keepdiag']), max_dev=1e-12)
