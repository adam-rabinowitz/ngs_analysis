'''concatFastq.py

Usage:
    
    concatFastq.py prefix <outfastq> <prefix>
    concatFastq.py specify <outfastq> <infastq>..
    
'''
# Import modules
import os
import sys
from general_python import moab, docopt, toolbox
from ngs_python.fastq import fastqFind
# Extract arguments
args = docopt.docopt(__doc__, version = 'v1')
# Find FASTQ files by prefix
if args['prefix']:
    indir, prefix = os.path.split(args['<prefix>'])
    print indir, prefix
    args['<infastq>'] = fastqFind.findFastq(prefix = prefix,
        dirList = [indir], pair = False)
    args['<infastq>'].sort()
# Check number of FASTQ files
if len(args['<infastq>']) < 2:
    sys.exit('\nCannot concatenate %s files\n' %(len(args['<infastq>'])))
# Check output file doesnt exist
if os.path.isfile(args['<outfastq>']):
    sys.exit('\nOutput file exists. No command submitted\n')
# Print input and out files
print '\nInput files:\n%s\n\nOutput file:\n%s\n' %(
    '\n'.join(args['<infastq>']), args['<outfastq>'])
# Get user response before concatenation
print "Enter 'concat' to concatenate: "
response = raw_input()
# Submit command
if response == 'concat':
    # Create and submit command
    command = 'zcat %s | gzip > %s' %(
        ' '.join(args['<infastq>']),
        args['<outfastq>']
    )
    moabID = moab.submitJob(command)
    print '\n%s\n' %(moabID)
else:
    print '\nNo command submitted\n'
