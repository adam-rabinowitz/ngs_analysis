'''meanCoverageCorrelation.py
    
Usage:
    meanCoverageCorrelation.py <intervals> <outprefix> <bam>... 
        [--minmap=<mm>] [--mincov=<mc>] [--rmdup] [--onebased]

Options:
    --minmap=<mm>  Minimum mapping quality for read [default: 0].
    --mincov=<mc>  Minimum coverage for interval [default: 5].
    --rmdup        Skip duplicate reads in calculating coverage.
    --onebased     Intervals have a one-based start. Otherwise a
        zero-based start is presumed.

'''
# Load required modules
import os
import numpy as np
from general_python import docopt, custom_plots
from ngs_python.bam import pysam_coverage
# Extract and check arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--minmap'] = int(args['--minmap'])
if args['--minmap'] < 0:
    raise ValueError('--minmap must be >= 0')
args['--mincov'] = int(args['--mincov'])
if args['--mincov'] < 1:
    raise ValueError('--mincov must be >= 1')
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
coverage = bamCov.mean_coverage(intervals=intervalList,
    map_quality=args['--minmap'], remove_dup=args['--rmdup'])
# Shorten bam names
for count, bam in enumerate(coverage):
    name = os.path.basename(bam)
    if name.endswith('.bam'):
        name = name[:-4]
    if count:
        nameList.append(name)
    else:
        nameList = [name]
coverage.columns = nameList
# Filter coverage and log convert
filtCov = coverage[coverage.min(axis=1) >= args['--mincov']]
log2Cov = np.log2(filtCov)
# Create output files
covFile = args['<outprefix>'] + '.mean_coverage.txt'
corFile = args['<outprefix>'] + '.correlation_matrix.txt'
corPlot = args['<outprefix>'] + '.correlation_matrix.png'
denPlot = args['<outprefix>'] + '.correlation_dendrogram.png'
# Save data to file
with open(corFile, 'w') as outfile:
    outfile.write('# BAM files:\n')
    for bam in args['<bam>']:
        outfile.write('#   {}\n'.format(bam))
    outfile.write('# Parameters:\n')
    outfile.write('#   min coverage - {}\n'.format(args['--mincov']))
    outfile.write('#   min map quality - {}\n'.format(args['--minmap']))
    outfile.write('#   remove duplicates - {}\n'.format(args['--rmdup']))
    outfile.write('#   supplied intervals - {}\n'.format(len(intervalList)))
    outfile.write('#   correlated intervals - {}\n'.format(len(filtCov.index)))
coverage.corr().to_csv(corFile, sep='\t', mode='append')
coverage.to_csv(covFile, sep='\t', index_label='interval')
custom_plots.correlationDendrogram(log2Cov, denPlot)
custom_plots.correlationMatrix(log2Cov, corPlot)
