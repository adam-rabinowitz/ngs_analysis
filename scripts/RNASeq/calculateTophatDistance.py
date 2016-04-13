''' calculateTophatDistance.py

Usage:
    
    calculateTophatDistance.py <nsortbam> [--maxdist=<maxdist>]
    
Options:
    
    --maxdist=<maxdist>  Maximum read fragment size [default: 2000]
    
'''
# Import required modules
import numpy as np
from general_python import docopt
from ngs_python.bam import pysamfunc
# Extract arguments and create list to store distances
args = docopt.docopt(__doc__,version = 'v1')
args['--maxdist'] = int(args['--maxdist'])
innerList = []
# Extract read pairs from BAM file
for read1, read2 in pysamfunc.pairGenerator(args['<nsortbam>']):
    # Extract fragment size
    outer, inner = pysamfunc.pairDist(
        read1 = read1,
        read2 = read2
    )
    # If outer is set reutn iner
    if outer and outer <= args['--maxdist']:
        innerList.append(inner)
# Calculate and print metrics
mean = np.mean(innerList, dtype=int)
stdv = np.std(innerList, dtype=int, ddof=1)
print mean, stdv
