import collections
import gzip

def paragraph_generator(
        path
    ):
    # Create output variable
    paragraph = []
    # Open input file and loop through file
    if path.endswith('.gz'):
        openfunc = gzip.open
    else:
        openfunc = open
    with openfunc(path) as infile:
        for line in infile:
            line = line.strip().decode('utf-8')
            # Append line to paragraph
            if line:
                paragraph.append(line)
            # Or return complete paragraphs
            else:
                if len(paragraph) > 0:
                    yield(paragraph)
                    paragraph = []
    # Yield final paragraph
    if len(paragraph) > 0:
        yield(paragraph)

def parse_terms(
        path
    ):
    # Create output variable and loop through go file
    output = collections.OrderedDict()
    for paragraph in paragraph_generator(path):
        # Check paragraphs of GO terms
        if paragraph[0] != '[Term]':
            continue
        if not paragraph[1].startswith('id: GO:'):
            raise ValueError(paragraph[1])
        if not paragraph[2].startswith('name: '):
            raise ValueError(paragraph[2])
        if not paragraph[3].startswith('namespace: '):
            raise ValueError(paragraph[3])
        # Extracrt data
        term = paragraph[1].split(' ', 1)[1]
        name = paragraph[2].split(' ', 1)[1]
        domain = paragraph[3].split(' ', 1)[1]
        output[term] = (name, domain)
    return(output)

terms = parse_terms('/g/furlong/genome/D.melanogaster/Dm6/6.13/go_annotation/go-basic.obo.gz')

def parse_genes(
        path, evidence=None
    ):
    # Check arguments
    if evidence is not None:
        if not isinstance(evidence, list):
            raise TypeError('evidence must be a list')
        for e in evidence:
            if not isinstance(e, str):
                raise TypeError('evidence must be a list of strings')
    # Create output variables
    genes2go = collections.defaultdict(set)
    go2genes = collections.defaultdict(set)
    if path.endswith('.gz'):
        openfunc = gzip.open
    else:
        openfunc = open
    with openfunc(path) as infile:
        for line in infile:
            # Decode line and skip headers
            line = line.decode('utf-8')
            if line.startswith('!'):
                continue
            # Extract data from line and check
            linedata = line.strip().split('\t')
            gene = linedata[1]
            if not gene.startswith('FBgn'):
                raise ValueError('Unexpected gene id: {}'.format(gene))
            term = linedata[4]
            if not term.startswith('GO:'):
                raise ValueError('Unexpected GO id: {}'.formar(term))
            # Extract and check evidence code
            code = linedata[6]
            if evidence is not None:
                if code not in evidence:
                    continue
            # Store data
            genes2go[gene].add(term)
            go2genes[term].add(gene)
    # Return data
    return(genes2go, go2genes)

def parse_results_file(
        path, sig, gene_header='gene', pvalue_header='padj'
    ):
    ''' Function to extract signifcant and insignificant genes from results
    file. Only returns genes with a reported pvalue.
    
    Args:
        path (str)- Path to results file.
        sig (float)- Threshold for gene significance.
        gene_header (str)- Column header for gene column of file.
        pvalue_header (str)- Column header for pvalue columns of file.
    
    Returns:
        allGenes - A set of all genes with a pvalue.
        sigGenes - A set of genes with a pvalue <= the supplied threshold.

    '''
    # Create output variables and loop through input file
    allGenes = set()
    sigGenes = set()
    with open(path, 'r') as infile:
        # Extract column for gene and pvalue from header
        header = infile.readline().strip().split('\t')
        pIndex = header.index('padj')
        gIndex = header.index('gene')
        # Find GO terms for genes
        for line in infile:
            # Extract gene name and pvalue or skip
            linedata = line.strip().split('\t')
            gene = linedata[gIndex]
            try:
                pvalue = float(linedata[pIndex])
            except ValueError:
                continue
            # Process genes with at least one associated go term
            allGenes.add(gene)
            if pvalue <= sig:
                sigGenes.add(gene)
    # Return data
    return(allGenes, sigGenes)

# Import required modules
import argparse
import gzip
import collections
import scipy.stats
import pandas
import statsmodels.stats.multitest as ssm

def go_analysis(
        allGenes, sigGenes, gene2GO, GO2gene, goAnno, minGO=5, maxGO=500,
        minSigGene=3
    ):
    ''' Function to perform gene ontology analaysis on results file from
    deseq. For genes to be included in the analysis they must be
    associated with a pvalue in the deseq results file and the must be
    associated with at least one gene ontology term in the gene2go file.
    
    Args:
        allgenes (set)- A set of all genes in the analysis.
        siggenes (set)- A set of all significant genes in the analysis.
        gene2go (dict)- A dictionary where the keys are gene IDs and the
            values are a set of gene ontology terms.
        go2genes (dict)- A dictionary where the keys are GO IDs and the
            values are a set of associated genes.
        mingo (int)- Minimum nuber of genes associated with a GO term.
        maxgo (int)- Maximum number of genes associated with a GO term.
        mingene (int)- Minimum number of significant genes associated with
            a GO term.
    
    '''
    # Check arguments
    print((minGO, maxGO, minSigGene))
    if not sigGenes.issubset(allGenes):
        raise ValueError('sigGenes must be a subset of allGenes')
    # Create variables to store set of genes associated with GO terms
    allGeneGO = set()
    sigGeneGO = set()
    allGeneGOdict = collections.defaultdict(set)
    sigGeneGOdict = collections.defaultdict(set)
    # Find genes with go terms
    for gene in allGenes:
        setGO = gene2GO[gene]
        if len(setGO) > 0:
            # Store data for genes
            allGeneGO.add(gene)
            for GO in setGO:
                allGeneGOdict[GO].add(gene)
            # Store data for significant genes
            if gene in sigGenes:
                sigGeneGO.add(gene)
                for GO in setGO:
                    sigGeneGOdict[GO].add(gene)
    # Set fixed values for hypegeometric mean
    N = len(allGeneGO)
    n = len(sigGeneGO)
    # Loop through GO terms for signifcant
    results = pandas.DataFrame()
    for GO, genes in sigGeneGOdict.items():
        # Set additional variables for hypergeometric test
        k = len(genes)
        K = len(allGeneGOdict[GO])
        # Skip genes with unwanted values
        if K > maxGO:
            continue
        if K < minGO:
            continue
        if k < minSigGene:
            continue
        # Calculate p-value
        pvalue = 1 - scipy.stats.hypergeom.cdf(k-1, N, K, n)
        goResults = pandas.Series([GO, goAnno[GO][0], goAnno[GO][1],
            K, N, k, n, pvalue])
        results[GO] = goResults
    # Maniuplate data frame and add column names
    results = results.transpose()
    results.columns = ['go_id', 'go_name', 'go_domain', 'K', 'N', 'k', 'n', 'pvalue']
    # Sort by pvalue, add fdr and return
    results.sort_values(by='pvalue', inplace=True)
    results['fdr'] = ssm.multipletests(results['pvalue'],method='fdr_bh')[1]
    return(results)

# Create variables to store gene and GO data
goAnno = parse_terms('/g/furlong/genome/D.melanogaster/Dm6/6.13/go_annotation/go-basic.obo.gz')
gene2GO, GO2gene = parse_genes('/g/furlong/genome/D.melanogaster/Dm6/6.13/go_annotation/gene_association.fb.gz')
allGenes, sigGenes = parse_results_file('/g/furlong/project/60_CTCF_RNASeq/data/deseq2/deseq2.results.txt', 0.05)
x = go_analysis(
    allGenes, sigGenes, gene2GO, GO2gene, goAnno, minGO=5, maxGO=500,
    minSigGene=3
)
x.to_csv('/g/furlong/project/60_CTCF_RNASeq/data/deseq2/deseq2.go_analysis.txt', sep='\t', index=False)





