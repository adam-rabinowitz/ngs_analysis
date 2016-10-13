"""translateResults.py

Usage:
    
    translateResults.py <results> <genecol> <translation> <outfile> 
        [--keepdup=<kkepdup>] [--noheader=<noheader>]
    
    translateResults.py (-h | --help)g
    
Options:
    
    --keepdup=<pa>          Keep genes with multiple translations.
    --noheader=<noheader>   Results file has no header.
    
"""
# Import required modules
from ngs_python.gtf import gene_conversion
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<genecol>'] = int(args['<genecol>'])
# Extract translations and generate output
tranDict, tranLog = gene_conversion.create_tran_dictionary(
    args['<translation>'])
gene_conversion.replace_gene_names(
    infile=args['<results>'], outfile=args['<outfile>'],
    genecol=args['<genecol>'], translation=tranDict,
    header=not(args['--noheader']), rmdup=not(args['--keepdup'])
)
