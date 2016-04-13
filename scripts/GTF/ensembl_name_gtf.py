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
nRE = re.compile('gene_name\s+"(.*?)";')
# Read in input file to gene dictionary
geneDict = collections.defaultdict(set)
with open(args['<gtffile>'], 'r') as inFile:
    for line in inFile:
        data = line.strip().split('\t')[8]
        ensembl = re.search(eRE, data).group(1)
        name = re.search(nRE, data). group(1)
        geneDict[ensembl].add(name)
# Create output file
counter = [0, 0, 0]
with open(args['<outfile>'], 'w') as outFile:
    for ensembl, names in geneDict.items():
        counter[0] += 1
        names = filter(None, names)
        if names:
            if len(names) > 1:
                counter[1] += 1
            for name in names:
                outFile.write('%s\t%s\n' %(ensembl, name))
        else:
            counter[2] += 1
            outFile.write('%s\t""\n' %(ensembl))
# Print counter
print '%s genes processed' %(counter[0])
print '%s genes with multiple names' %(counter[1])
print '%s genes with no names' %(counter[2])
