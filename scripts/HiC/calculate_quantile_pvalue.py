'''calc_quantile_pvalue.py

Usage:
    
    calc_quantile_pvalue.py <samples1> <samples2> <names> <quantiles>
        <indir> <outfile>
    
Options:
    -h --help  Show this screen.

'''
import os
from general_python import docopt
from ngs_python.structure import analyseInteraction
args = docopt.docopt(__doc__,version = 'v1')
# Split comma sperated arguments
args['<samples1>'] = args['<samples1>'].split(',')
args['<samples2>'] = args['<samples2>'].split(',')
args['<names>'] = args['<names>'].split(',')
args['<quantiles>'] = args['<quantiles>'].split(',')
args['<quantiles>'] = [float(x) for x in args['<quantiles>']]
# Check input and output directories exist
if not os.path.isdir(args['<indir>']):
    raise IOError('{} not found'.format(args['<indir>']))
if not os.path.isdir(os.path.dirname(args['<outfile>'])):
    raise IOError('{} not found'.format(os.path.dirname(args['<outfile>'])))
# Extract matrix dictionary
analyseMatrices = analyseInteraction.compare_paired_matrices(
    samples1=args['<samples1>'], samples2=args['<samples2>'],
    indir=args['<indir>'], conditions=args['<names>'],
    suffix='normMatrix.gz'
)
# Log matrices in output file
matrices = analyseMatrices.matrices
with open(args['<outfile>'], 'w') as outfile:
    for cond in matrices:
        outfile.write('# {} prefixes:\n'.format(cond))
        for smpl in matrices[cond]:
            outfile.write('#    {}\n'.format(smpl))
    for cond in matrices:
        outfile.write('# {} matrices:\n'.format(cond))
        for smpl in matrices[cond]:
            for matrix in matrices[cond][smpl]:
                outfile.write('#    {}\n'.format(matrix))
# Calculate distance p-values
quanPvalue = analyseMatrices.calculate_quantile_pvalue(args['<quantiles>'])
# Print output dataframe
quanPvalue.to_csv(args['<outfile>'], sep='\t', mode = 'a',
    index_label='quantile')
