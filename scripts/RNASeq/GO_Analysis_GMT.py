"""GO_Analysis_GMT.py

Usage:
    
    GO_Analysis_GMT.py <genefile> <gmtfile> <outfile> [--mingo=<mingo>]
        [--maxgo=<maxgo>] [--mingene=<mingene>]
    
    GO_Analysis_GMT.py (-h | --help)
    
Options:
    
    --mingo=<mingo>      Minimum number of genes for GO term [default: 5]
    --maxgo=<maxgo>      Maximum number of genes for GO term [default: 1000]
    --mingene=<mingene>  Minimum number of significant genes [default: 5]
    --help               Output this message
    
"""
# Import required modules
import argparse
import gzip
import collections
import scipy.stats
import pandas as pd
import statsmodels.stats.multitest as ssm
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--mingo'] = int(args['--mingo'])
args['--maxgo'] = int(args['--maxgo'])
args['--mingene'] = int(args['--mingene'])
# Create variables to store gene and GO data
genes2GO = collections.defaultdict(set)
GO2genes = collections.defaultdict(set)
GO2anno = {}
# Create variables to store data from differential expression
allGeneSet = set()
sigGeneSet = set()
sigGenesGO = collections.defaultdict(set)
# Loop through GO data and populate gene-GO dictionaries
with open(args['<gmtfile>'], "r") as ifile:
    # Associate GO terms with gene names
    for line in ifile:
        # Extract data
        lineData = line.strip().split('\t')
        GO, anno = lineData[:2]
        GO2anno[GO] = anno
        # Extract GO terms associated with genes
        for gene in lineData[2:]:
            allGeneSet.add(gene)
            genes2GO[gene].add(GO)
            GO2genes[GO].add(gene)
# Loop through differential expression and find GO terms
with open(args['<genefile>'], 'r') as ifile:
    # Find GO terms for genes
    for line in ifile:
        # Extract gene name and associate go terms
        gene = line.strip().split('\t')[0]
        goSet = genes2GO[gene]
        # Process genes with associated GO terms
        if goSet:
            # Add gene to list of annotated genes
            sigGeneSet.add(gene)
            for GO in goSet:
                sigGenesGO[GO].add(gene)

###############################################################################
## Process and save results
###############################################################################
results = pd.DataFrame()
# Set fixed values for hypegeometric mean
N = len(allGeneSet)
n = len(sigGeneSet)
# Loop through GO terms for signifcant
for GO, genes in sigGenesGO.iteritems():
    # Set additional variables for hypergeometric test
    k = len(genes)
    K = len(GO2genes[GO])
    if (K >= args['--mingo'] and
        K <= args['--maxgo'] and
        k >= args['--mingene']):
        # Calculate p-value
        pvalue = 1 - scipy.stats.hypergeom.cdf(
            k-1,N,K,n
        )
    else:
        pvalue = pd.np.nan
    goResults = pd.Series([GO, GO2anno[GO], K, N, k, n, pvalue])
    results[GO] = goResults
# Maniuplate data frame and add column names
results = results.transpose()
results.columns = ['go_term', 'go_anno', 'K', 'N', 'k', 'n', 'pvalue']
# Sort by pvalue and add fdr
results.sort_values(by='pvalue', inplace=True)
pvalueIndex = pd.np.where(~pd.isnull(results['pvalue']))
pvalue = results['pvalue'].iloc[pvalueIndex]
fdr = ssm.multipletests(pvalue, method='fdr_bh')[1]
fdrSeries = pd.Series(index = results.index)
fdrSeries.iloc[pvalueIndex] = fdr
results['fdr'] = fdrSeries
# Save results to file
results.to_csv(args['<outfile>'], sep = '\t', index = False)
