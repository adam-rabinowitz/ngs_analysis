'''filterGTF.py

Usage:
    
    filterGTF.py <ingtf> <fastafile> <outgtf>

'''
# Import required modules
from general_python import docopt
import re
import os
# Extract arguments
args = docopt.docopt(__doc__, version='1.0')
if not os.path.isfile(args['<ingtf>']):
    raise IOError('Input GTF file not found')
if not os.path.isfile(args['<fastafile>']):
    raise IOError('Input FASTA file not found')
if os.path.isfile(args['<outgtf>']):
    raise IOError('Output GTF file already exists')
# Extract chromosome names from fasta file
chromList = []
with open(args['<fastafile>']) as fastafile:
    for line in fastafile:
        if line.startswith('>'):
            chrom = re.search('^>\s*([^\s]+)', line).group(1)
            chromList.append(chrom)
print 'Chromosomes:\n\t%s' %('\n\t'.join(chromList))
# Extract arguments
args = docopt.docopt(__doc__, version='1.0')
# Create output gtf file
count = [0,0,0,0]
with open(args['<ingtf>']) as ingtf:
    with open(args['<outgtf>'], 'w') as outgtf:
        for line in ingtf:
            # Write header lines
            if line.startswith('#'):
                outgtf.write(line)
                continue
            # Split data lines and count total
            gtfdata = line.split('\t')
            count[0] += 1
            # Find acceptable entries
            if gtfdata[0] in chromList and gtfdata[2] == 'exon':
                outgtf.write(line)
                count[3] += 1
            # Count entries on unwanted chromosome
            elif gtfdata[2] == 'exon':
                count[1] += 1
            # Count non-exonic entries
            else:
                count[2] += 1
print 'Metrics:'
print '\tTotal: %s' %(count[0])
print '\tUnwanted Chromosome: %s' %(count[1])
print '\tNot An Exon: %s' %(count[2])
print '\tAccepted: %s' %(count[3])
