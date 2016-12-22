'''coverageHistogram.py
    
Usage:
    coverageHistogram.py <intervals> <outfile> <bam>... 
        [--minmap=<mm>] [--maxcov=<mc>] [--rmdup] [--rmsec] [--onebased]

Options:
    --minmap=<mm>  Minimum mapping quality for read [default: 0].
    --maxcov=<mc>  Maximum coverage to calculate [default: 500].
    --rmdup        Skip duplicate reads in calculating coverage.
    --rmsec        Skip secondary reads in calculating coverage.
    --onebased     Intervals have a one-based start. Otherwise a
        zero-based start is presumed.

'''
# Load required modules
from general_python import docopt
from ngs_python.bam import pysam_coverage
import pandas as pd
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--minmap'] = int(args['--minmap'])
args['--maxcov'] = int(args['--maxcov'])
# Open interval list file and extract data
intervalList = []
with open(args['<intervals>']) as intervalFile:
    for line in intervalFile:
        chrom, start, end = line.strip().split('\t')[:3]
        intervalList.append((chrom, int(start), int(end)))
# Adjust intervals if they are one based
if args['--onebased']:
    intervalList = [(x[0], x[1] - 1, x[2]) for x in intervalList]
# Create output dataframe
outDF = pd.DataFrame(
    index=range(args['--maxcov'] + 1),
    columns=args['<bam>'])
# Extract mean coverage for intervals
for bam in args['<bam>']:
    covCalc = pysam_coverage.single_coverage(bam)
    outDF[bam] = covCalc.coverage_histogram(
        intervals=intervalList, max_cov=args['--maxcov'],
        map_quality=args['--minmap'], remove_dup=args['--rmdup'],
        remove_secondary=args['--rmsec'])
# Save intervals to file
if args['<outfile>'].endswith('.gz'):
    compression = 'gzip'
else:
    compression = None
outDF.to_csv(args['<outfile>'], sep='\t', index_label='coverage',
    compression=compression)
