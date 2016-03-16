"""GO_Analysis.py

Usage:
    
    GO_Analysis.py <diffile> <gofile> <outfile> [--mingo=<mingo>]
        [--maxgo=<maxgo>] [--mingene=<mingene>] [--padj=<padj>]
    
    GO_Analysis.py (-h | --help)g
    
Options:
    
    --mingo=<mingo>      Minimum number of genes for GO term [default: 5]
    --maxgo=<maxgo>      Maximum number of genes for GO term [default: 1000]
    --mingene=<mingene>  Minimum number of significant genes [default: 5]
    --padj=<padj>        Adjusted pvalue for significant genes [default: 0.05]
    --help               Output this message
    
"""
# Import required modules
import argparse
import gzip
import collections
import scipy.stats
import pandas
import statsmodels.stats.multitest as ssm
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--mingo'] = int(args['--mingo'])
args['--maxgo'] = int(args['--maxgo'])
args['--mingene'] = int(args['--mingene'])
args['--padj'] = float(args['--padj'])
# Create variables to store gene and GO data
goAnnotation = {}
genes2GO = collections.defaultdict(set)
GO2genes = collections.defaultdict(set)
# Create variables to store data from differential expression
geneSet = set()
sigGeneSet = set()
allGenesGO = collections.defaultdict(set)
sigGenesGO = collections.defaultdict(set)
# Loop through GO data and populate gene-GO dictionaries
with gzip.open(args['<gofile>'], "r") as ifile:
    # Skip header
    next(ifile)
    # Associate GO terms with gene names
    for line in ifile:
        # Extract date
        gene, GO, name, evidence, domain = line.strip('\n').split('\t')
        # Extract GO terms associated with genes
        genes2GO[gene].add(GO)
        # Extract genes associated with GO terms
        GO2genes[GO].add(gene)
        # Extract names associated with GO terms
        goAnnotation[GO] = (name, domain)
# Loop through differential expression and find GO terms
with open(args['<diffile>'], 'r') as ifile:
    # Extract column for gene and adjusted pvalue from header
    header = ifile.readline().strip().split('\t')
    padjIndex = header.index('padj')
    # Find GO terms for genes
    for line in ifile:
        # Extract gene name and adjuste p-value from data
        data = line.strip().split('\t')
        gene = data[0]
        try:
            padj = float(data[padjIndex])
        except ValueError:
            padj = 1.0
        # Extract GO terms for gene
        goList = genes2GO[gene]
        # Process genes with associated GO terms
        if goList:
            # Add gene to list of annotated genes
            geneSet.add(gene)
            # Store gene names associated with GO terms
            if padj <= args['--padj']:
                # Add gene to list of significant genes
                sigGeneSet.add(gene)
                for GO in goList:
                    sigGenesGO[GO].add(gene)
                    allGenesGO[GO].add(gene)
            else:
                for GO in goList:
                    allGenesGO[GO].add(gene)

################################################################################
## Process and save results
################################################################################
results = pandas.DataFrame()
# Set fixed values for hypegeometric mean
N = len(geneSet)
n = len(sigGeneSet)
# Loop through GO terms for signifcant
for GO, genes in sigGenesGO.iteritems():
    # Set additional variables for hypergeometric test
    k = len(genes)
    K = len(allGenesGO[GO])
    if (K >= args['--mingo'] and
        K <= args['--maxgo'] and
        k >= args['--mingene']):
        # Calculate p-value
        pvalue = 1 - scipy.stats.hypergeom.cdf(
            k-1,N,K,n
        )
        goResults = pandas.Series([GO, goAnnotation[GO][0], goAnnotation[GO][1],
            K, N, k, n, pvalue])
        results[GO] = goResults
# Maniuplate data frame and add column names
results = results.transpose()
results.columns = ['go_term', 'go_name', 'go_domain', 'K', 'N', 'k', 'n', 
    'pvalue']
# Sort by pvalue and add fdr
results.sort_values(by='pvalue', inplace=True)
results['fdr'] = ssm.multipletests(results['pvalue'],method='fdr_bh')[1]
# Save results to file
results.to_csv(args['<outfile>'], sep = '\t', index = False)
