'''binOverlapCount.py
    
Usage:
    binOverlapCount.py <bam> <binsize> <overlap> <mapq> <outfile>
        [--equal] [--rmdup] [--rmsec] [--rmsup]

Options:
    --equal  Bins should be equally sized.
    --rmdup  Remove duplicate reads.
    --rmsec  Remove secondary reads.
    --rmsup  Remove supplementray reads.

'''
# Load required modules
import os
from ngs_python.bam import pysam_coverage
from general_python import docopt
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<mapq>'] = int(args['<mapq>'])
args['<binsize>'] = int(args['<binsize>'])
args['<bam>'] = os.path.abspath(args['<bam>'])
args['<outfile>'] = os.path.abspath(args['<outfile>'])
# Create coverage object and extract data
sc = pysam_coverage.single_coverage(args['<bam>'])
binDict = sc.count_bin_overlaps(
    binSize=int(args['<binsize>']), binEqual=args['--equal'],
    mapQ=args['<mapq>'], overlap=args['<overlap>'], rmDup=args['--rmdup'],
    rmSec=args['--rmsec'], rmSup=args['--rmsup']
)
# Open outfile and print parameters
outfile = open(args['<outfile>'], 'w')
outfile.write('# input file: {}\n'.format(args['<bam>']))
outfile.write('# bin size: {}\n'.format(args['<binsize>']))
outfile.write('# bin size equal: {}\n'.format(args['--equal']))
outfile.write('# ovarlap type: {}\n'.format(args['<overlap>']))
outfile.write('# minimum map quality: {}\n'.format(args['<mapq>']))
outfile.write('# remove duplicate alignments: {}\n'.format(args['--rmdup']))
outfile.write('# remove secondary alignments: {}\n'.format(args['--rmsec']))
outfile.write('# remove supplementary alignments: {}\n'.format(args['--rmsup']))
# Add output data
outfile.write('chr\tstart\tend\tcount\n')
for chrom in binDict:
    chrData = binDict[chrom]
    for bindata in zip(
                chrData['start'], chrData['end'], chrData['count']
            ):
            output = '{}\t{}\t{}\t{}\n'.format(
                chrom, bindata[0] + 1, bindata[1], bindata[2])
            outfile.write(output)
outfile.close()
