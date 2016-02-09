'''alignSortIndex.py

Usage:

    alignSortIndex.py <sampledata> <indir> <outdir> <index>
        [--threads=<threads>] [--illumina] [--paired]
    
Options:

    --threads=<threads>  Number of threads [default: 1]
    --illumina           Illumina format file names
    --pair=<pair>        Input reads should be paired
    --help               Output this message

'''
# Import required modules
from general_functions import docopt
from ngs_python.fastq import fastqFind
print 'hello'
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
print args
# Check numerical arguments
args['--threads'] = int(args['--threads'])
# Split sample data
args['name'], args['prefix'] = args['<sampledata>'].split(',')
print args
