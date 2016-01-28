"""HiC_BinCounts.py

Usage:
    
    concatIlluminaFrags.py <indir> <outdir> [--unpaired]
    
    concatIlluminaFrags.py (-h | --help)
    
Options:
    
    --unpaired  Do not check for paired reads
    
"""

# Import required modules
import os
import re
import collections
from general_functions import docopt
from general_functions import moab
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Create regular expression to find fastq files
rePattern = re.compile('(.*?)_L\\d{3}_(R[12])_\\d{3}\\.fastq\\.gz$')
# Create dictionary
fastqDict = {}
# Extract input files
for infile in os.listdir(args['<indir>']):
    # Match fastq files with regular expression
    match = re.match(rePattern, infile)
    if not match:
        continue
    # Add files to output dictionary
    name = match.group(1)
    read = match.group(2)
    if name not in fastqDict:
        fastqDict[name] = {'R1':[],'R2':[]}
    fastqDict[name][read].append(infile)
# Process samples
for sample in fastqDict:
    # Sort and extract reads
    fastqDict[sample]['R1'].sort()
    R1 = fastqDict[sample]['R1']
    fastqDict[sample]['R2'].sort()
    R2 = fastqDict[sample]['R2']
    # Check for pairs
    if not args['--unpaired']:
        print 'Yo'
        # Check equal number of read1 and read2 files
        if len(R1) != len(R2):
            raise IOError('Unpaired reads found for %s' %(sample))
        # Check for matching read1 and read2 files
        for number, r1file in enumerate(R1):
            r2file = re.sub('R1(?=_\\d{3}\\.fastq\\.gz$)', 'R2', r1file)
            if R2[number] != r2file:
                raise IOError('Unpaired files found for %s %s' %(r1file,
                    r2file))
    # Create command
    for infiles in [R1, R2]:
        # Create output file
        outfile = re.sub('_\\d{3}(?=\\.fastq\\.gz$)', '', infiles[0])
        outfile = args['<outdir>'] + outfile
        # Append path to input files
        infiles = [args['<indir>'] + r for r in infiles]
        # Create command
        command = 'zcat %s | gzip > %s' %(' '.join(infiles), outfile)
        moabID = moab.submitjob(command, processor=2)
        print moabID
        
