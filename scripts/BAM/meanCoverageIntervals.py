'''meanCoverageIntervals.py
    
Usage:
    meanCoverageIntervals.py <intervals> <outfile> <bam>... 
        [--minmap=<mm>] [--skipdup] [--onebased]

Options:
    --minmap=<mm>  Minimum mapping quality for read [default: 0].
    --skipdup      Skip duplicate reads in calculating coverage.
    --onebased     Intervals have a one-based index. Otherwise a zero based
        index is presumed.

'''
# Load required modules
import pandas as pd
import pysam
from general_python import docopt
from ngs_python.bam import pysam_coverage
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--minmap'] = int(args['--minmap'])
# Check all BAM files are aligned to the same reference genome
for count, bam in enumerate(args['<bam>']):
    # Extract chromsome names and lengths
    bamFile = pysam.AlignmentFile(bam)
    sampleDict = dict(zip(bamFile.references, bamFile.lengths))
    bamFile.close()
    # Compare lengths and names between samples
    if count:
        if sampleDict != refDict:
            raise ValueError('Not all BAM files have the same reference')
    else:
        refDict = sampleDict
# Create variables to store interval data
intervalList = []
nameList = []
# Open interval list file and extract data
with open(args['<intervals>']) as intervalFile:
    for line in intervalFile:
        chrom, start, end = line.strip().split('\t')[:3]
        # Create and store interval name
        name = '{}:{}-{}'.format(chrom, start, end)
        nameList.append(name)
        # Create and store intervals
        if args['--onebased']:
            start -= 1
        intervalList.append((chrom, int(start), int(end)))
# Create output dataframe
outDF = pd.DataFrame(columns=args['<bam>'], index=nameList)
for count, bam in enumerate(args['<bam>']):
    count = pysam_coverage.coverage(bam)
    if count == 0:
        count.check_intervals(intervalList)
    cov = count.mean_coverage(
        intervals = intervalList,
        map_quality = args['--minmap'],
        remove_dup = args['--skipdup'],
        check_intervals = False
    )
    outDF[bam] = cov
outDF.to_csv(args['<outfile>'], sep='\t')
