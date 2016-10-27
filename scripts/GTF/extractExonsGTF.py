'''extractExonsGTF.py
    
Usage:
    
    extractExonsGTF.py (gene|tran|name) <id> <gtf> <bed>

'''
from ngs_python.gtf import extract_gtf
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Create output dataframe
if args['gene']:
    ident_type = 'gene'
elif args['tran']:
    ident_type = 'tran'
elif args['name']:
    ident_type = 'name'
# Create bed file
extract_gtf.extract_exons_bed(
    args['<id>'],
    ident_type,
    args['<gtf>'],
    args['<bed>']
)
