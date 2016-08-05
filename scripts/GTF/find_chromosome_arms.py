'''find_chromsome_arms.py

Usage:
    
    find_chromosome_arms.py <fasta> <gff> <outfile> [--bed] [--name]

'''
# Import required modules
import collections
import re
from general_python import docopt
from Bio import SeqIO
# Extract arguments
args = docopt.docopt(__doc__, version='1.0')
# Store lengths of chromosomes
lengths = collections.OrderedDict()
with open(args['<fasta>'], "r") as fasta:
    for record in SeqIO.parse(fasta, "fasta"):
        lengths[record.id] = len(record.seq)
# Extract centromeres
centro = {}
with open(args['<gff>'], 'r') as gff:
    count = 0
    for line in gff:
        count += 1
        if line.startswith('#'):
            continue
        lineData = line.strip().split('\t')
        try:
            two = lineData[2]
            eight = lineData[8]
        except IndexError:
            continue
        if (two == 'centromere'
            and 'ID=CEN' in eight):
            chrom = lineData[0][3:]
            start = int(lineData[3])
            end = int(lineData[4])
            if chrom in centro:
                raise ValueError('Two {} centromeres identified'.format(chrom))
            else:
                centro[chrom] = (start, end)
print(centro)
print(lengths)
# Extract intervals
with open(args['<outfile>'], 'w') as outfile:
    for chrom in lengths:
        if chrom in centro:
            rEnd = lengths[chrom]
            lEnd = centro[chrom][0] - 1
            if args['--bed']:
                lStart = 0
                rStart = centro[chrom][1]
            else:
                lStart = 1
                rStart = centro[chrom][1] + 1 
            lArm = [chrom, lStart, lEnd]
            rArm = [chrom, rStart, rEnd]
            if args['--name']:
                lArm.append(chrom + '_ArmL')
                rArm.append(chrom + '_ArmR')
            output = '{}\n{}\n'.format(
                '\t'.join(map(str, lArm)),
                '\t'.join(map(str, rArm))
            )
            outfile.write(output)
        else:
            print 'No Centromere for chromsosome {}'.format(chrom)
