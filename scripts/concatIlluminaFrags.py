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
from general_python import docopt, moab
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Create regular expression to find fastq files
rePattern = re.compile('(.*?_L\\d{3}_R[12])_\\d{3}\\.fastq\\.gz$')
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
    path = os.path.join(args['<indir>'], infile)
    if name in fastqDict:
        fastqDict[name].append(path)
    else:
        fastqDict[name] = [path]
# Process samples
for name in fastqDict.keys():
    # Process unpaired reads
    if args['--unpaired']:
        # Extract and sort reads
        reads = fastqDict.pop(name)
        reads.sort()
        # Create outfile and command
        outfile = os.path.join(args['<outdir>'], name + '.fastq.gz')
        command = 'zcat %s | gzip > %s' %(' '.join(reads), outfile)
        # Check for missing files
        for count, rd in enumerate(reads):
            # Check for missing files
            number = int(rd[-12:][:3])
            if (number - count) != 1:
                raise IOError('Missing fragments found for %s' %(name))
    # Process paired reads
    else:
        # Extract files
        if name.endswith('1'):
            read1 = fastqDict.pop(name)
            read2 = fastqDict.pop(name[:-1] + '2')
        else:
            continue
        # Sort read files
        read1.sort()
        read2.sort()
        # Create outfiles
        outfile1 = os.path.join(args['<outdir>'], name[:-1] + '1.fastq.gz')
        outfile2 = os.path.join(args['<outdir>'], name[:-1] + '2.fastq.gz')
        # Check equal number of read1 and read2 files
        if len(read1) != len(read2):
            raise IOError('Unpaired reads found for %s' %(name))
        for count, (R1, R2) in enumerate(zip(read1, read2)):
            # Check pairing
            if  R1[-12:] != R1[-12:]:
                raise IOError('Unpaired fragments found for %s' %(name))
            # Check for missing files
            number = int(R1[-12:][:3])
            if (number - count) != 1:
                raise IOError('Missing fragments found for %s' %(name))
        # Create commands
        command1 = 'zcat %s | gzip > %s' %(' '.join(read1), outfile1)
        command2 = 'zcat %s | gzip > %s' %(' '.join(read2), outfile2)
        moab.submitJob(command1)
        moab.submitJob(command2)
        print '%s paired files concatenated to prefix %s' %(len(read1), outfile1[:-10])
        
