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
gRE = re.compile('gene_id\s+"(.*?)";')
tRE = re.compile('transcript_id\s+"(.*?)";')
# Read in input file to gene dictionary
geneDict = collections.defaultdict(set)
with open(args['<gtffile>'], 'r') as inFile:
    for line in inFile:
        data = line.strip().split('\t')[8]
        gene = re.search(gRE, data).group(1)
        tran = re.search(tRE, data).group(1)
        geneDict[gene].add(tran)
# Create output file
counter = [0, 0, 0]
with open(args['<outfile>'], 'w') as outFile:
    for gene, transcripts in geneDict.items():
        counter[0] += 1
        transcripts = filter(None, transcripts)
        if transcripts:
            if len(transcripts) > 1:
                counter[1] += 1
            for tran in transcripts:
                outFile.write('%s\t%s\n' %(gene, tran))
        else:
            counter[2] += 1
            outFile.write('%s\t""\n' %(gene))
# Print counter
print '%s genes processed' %(counter[0])
print '%s genes with multiple transcripts' %(counter[1])
print '%s genes with no transcript' %(counter[2])
