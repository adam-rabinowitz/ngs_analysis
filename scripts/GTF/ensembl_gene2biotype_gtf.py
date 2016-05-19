'''ensembl_name_gtf.py

Usage:
    
    ensembl_name_gtf.py <gtffile> <outfile>

'''
# Import required modules
import collections
import re
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__, version='1.0')
# Create regular expressions
eRE = re.compile('gene_id\s+"(.*?)";')
bRE = re.compile('gene_biotype\s+"(.*?)";')
# Read in input file to gene dictionary
geneDict = collections.defaultdict(set)
with open(args['<gtffile>'], 'r') as inFile:
    for line in inFile:
        if line.startswith('#'):
            continue
        data = line.strip().split('\t')[8]
        ensembl = re.search(eRE, data).group(1)
        biotype = re.search(bRE, data).group(1)
        geneDict[ensembl].add(biotype)
# Create output file
counter = [0, 0, 0]
with open(args['<outfile>'], 'w') as outFile:
    for ensembl, biotype in geneDict.items():
        counter[0] += 1
        biotype = filter(None, biotype)
        if biotype:
            if len(biotype) > 1:
                counter[1] += 1
            for bio in biotype:
                outFile.write('%s\t%s\n' %(ensembl, bio))
        else:
            counter[2] += 1
            outFile.write('%s\t""\n' %(ensembl))
# Print counter
print '%s genes processed' %(counter[0])
print '%s genes with multiple biotype' %(counter[1])
print '%s genes with no biotype' %(counter[2])
