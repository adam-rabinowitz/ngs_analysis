"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <regions> <mincount> <outdir> <infiles>...
        [--path=<paths>] [--threads=<threads>] [--memory=<memory>]
        [--iter=<iter>] [--keepdiag]
    
    generateNormalisedMatrix.py (-h | --help)
    
Options:
    
    --path=<path>        Path to ic_mes executable [default: ic_mes]
    --threads=<threads>  Number of threads [default: 1]
    --memory=<memory>    Memory (MB) per thread [default: 2000]
    --iter=<iter>        Iterations of ic_mes [default: 1000]
    --keepdiag           Keep diagonal of matrix, othersise set to zero.
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
    matrixList=args['<infiles>'], regionFile=args['<regions>'],
    rmDiag=not(args['--keepdiag']))
# Normalise sub matrices
normMatrices.normaliseSubMatrices(
    outDir=args['<outdir>'], minCount=args['<mincount>'], path=args['--path'],
    threads=args['--threads'], memory=args['--memory'],
    iterations=args['--iter'])
