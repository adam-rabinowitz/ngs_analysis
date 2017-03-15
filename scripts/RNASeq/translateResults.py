"""translateResults.py

Usage:
    
    translateResults.py <results> <genecol> <translation> <outfile> 
        [--keepdup] [--noheader]
    
    translateResults.py (-h | --help)
    
Options:
    
    --keepdup    Keep genes with multiple translations.
    --noheader   Results file has no header.
    
"""
# Import required modules
from ngs_python.gtf import gene_conversion
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<genecol>'] = int(args['<genecol>'])
# Extract translations and print log
tranDict, tranLog1 = gene_conversion.create_tran_dictionary(
    args['<translation>'])
tranLog2, tranLog2 = gene_conversion.replace_gene_names(
    infile=args['<results>'], outfile=args['<outfile>'],
    genecol=args['<genecol>'], translation=tranDict,
    header=not(args['--noheader']), rmdup=not(args['--keepdup'])
)
