'''meanCoverageIntervals.py
    
Usage:
    meanCoverageIntervals.py <intervals> <outfile> <bam>... 
        [--minmap=<mm>] [--rmdup] [--onebased] [--header]

Options:
    --minmap=<mm>  Minimum mapping quality for read [default: 0].
    --rmdup        Skip duplicate reads in calculating coverage.
    --onebased     Intervals have a one-based start. Otherwise a
        zero-based start is presumed.

'''
# Load required modules
from general_python import docopt
from ngs_python.bam import pysam_coverage
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--minmap'] = int(args['--minmap'])
# Open interval list file and extract data
intervalList = []
with open(args['<intervals>']) as intervalFile:
    for line in intervalFile:
        chrom, start, end = line.strip().split('\t')[:3]
        intervalList.append((chrom, int(start), int(end)))
# Adjust intervals if they are one based
if args['--onebased']:
    intervalList = [(x[0], x[1] - 1, x[2]) for x in intervalList]
# Extract mean coverage for intervals
bamCov = pysam_coverage.multiple_coverage(args['<bam>'])
outDF = bamCov.mean_coverage(intervals=intervalList,
    map_quality=args['--minmap'], remove_dup=args['--rmdup'])
# Save intervals to file
if args['<outfile>'].endswith('.gz'):
    compression = 'gzip'
else:
    compression = None
outDF.to_csv(args['<outfile>'], sep='\t', index_label='interval',
    compression=compression)
