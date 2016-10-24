import collections
import re
import pandas as pd
import numpy as np
from scipy.stats import hypergeom as ssh
import statsmodels.stats.multitest as ssm

def create_tran_dictionary(tranfile):
    ''' Parses text files to create dictionary translating
    between two sets of names.
    
    Args:
        tranfile - Full path to tab delimited text file listing input
        names in the first column and output names in second column.
    
    Returns:
        tranDict - A dictionary where keys are input names and values
        are lists of output names.
        log - A dictionary listing the number of output names found
        for the input names.
    
    '''
    # Sequentially process lines in file and add to dictionary
    tranDict = {}
    with open(tranfile) as infile:
        for line in infile:
            # Extract input and output genes from line
            lineData = line.strip().split('\t')
            if len(lineData) == 1:
                inGene = lineData[0]
                outGene = ''
            elif len(lineData) == 2:
                inGene, outGene = lineData
            else:
                raise ValueError('Line does not have 1 or 2 elements')
            # Add gene data to dictionary
            if inGene in tranDict:
                if outGene:
                    tranDict[inGene].add(outGene)
            else:
                tranDict[inGene] = set()
                if outGene:
                    tranDict[inGene].add(outGene)
    # Convert sets to lists
    tranDict = {k:list(v) for k,v in tranDict.iteritems()}
    # Create log and return data
    log = collections.defaultdict(int)
    for values in tranDict.itervalues():
        log[len(values)] += 1
    return(tranDict, log)

def replace_gene_names(
        infile, outfile, genecol, translation, header, rmdup=True
    ):
    ''' Parses input text files to create output text file where input
    gene names are replaced by gene names in translation dictionary.
    
    Args:
        infile (str)- Full path to input file.
        outfile (str)- Full path to output file.
        genecol (int)- Column in input file containing gene names.
        translation (dict)- Dictionary containing gene translations.
        header (bool)- Does input file contain a header.
        rmdup (bool)- Remove genes with multiple translations.
    
    Returns:
        log - A dictionary listing number of genes
    
    '''
    # Check arguments
    if not isinstance(genecol, int) or genecol < 0:
        raise ValueError('Unacceptable value for geneCol')
    if not isinstance(translation, dict):
        raise ValueError('Unacceptable value for translation')
    if not isinstance(header, bool):
        raise ValueError('Unacceptable value for header')
    if not isinstance(rmdup, bool):
        raise ValueError('Unacceptable value for rmdup')
    # Check names
    outGeneCounter = collections.defaultdict(int)
    with open(infile, 'r') as filein:
        if header:
            filein.next()
        for line in filein:
            lineData = line.strip().split('\t')
            inGene = lineData[genecol]
            outGenes = translation[inGene]
            for outGene in outGenes:
                outGeneCounter[outGene] += 1
    # Create output files
    log = {'total':0, 'multiOut':0, 'multiIn':0, 'noOut':0, 'single':0}
    with open(infile, 'r') as filein, open(outfile, 'w') as fileout:
        # Write header
        if header:
            line = filein.next()
            fileout.write(line)
        # Loop through line and create output
        for line in filein:
            # Parse line and extract input and output genes
            log['total'] += 1
            lineData = line.strip().split('\t')
            inGene = lineData[genecol]
            outGenes = translation[inGene]
            # Count and skip input genes with no output
            if len(outGenes) == 0:
                log['noOut'] += 1
                continue
            # Count and skip input genes with multiple output
            if rmdup and len(outGenes) > 1:
                log['multiOut'] += 1
                continue
            # Loop through output genes
            for outGene in outGenes:
                # Count and skip output genes with multiple input
                if rmdup and outGeneCounter[outGene] > 1:
                    log['multiIn'] += 1
                    continue
                # Write output file
                log['single'] += 1
                lineData[genecol] = outGene
                fileout.write('{}\n'.format('\t'.join(lineData)))
    # Return log
    return(log)

def extract_gene_results(
        results, geneCol, statCol, statMax, header
    ):
    ''' Parses text files to extract significant genes for ontology
    analysis.
    
    Args:
         results - Full path to text file.
         geneCol (int) - Column containing gene names.
         statCol (int) - Column containing significance statistic.
         statMax (float) - Threshold of significance for statistic.
         header (bool) - Whether text file contains a header.
    
    Returns:
        allGenes (set)- The set of all genes in file with an
            associated statistic.
        sigGenes (set)- The set of genes with a significant statistic.
    
    '''
    # Check arguments
    if not isinstance(geneCol, int) or geneCol < 0:
        raise ValueError('Unacceptable value for geneCol')
    if not isinstance(statCol, int) or statCol < 0:
        raise ValueError('Unacceptable value for statCol')
    if not isinstance(statMax, float) or 0 > statMax > 1:
        raise ValueError('Unacceptable value for statMax')
    if not isinstance(header, bool):
        raise ValueError('Unacceptable value for header')
    # Create variables to store data
    allGenes = set()
    sigGenes = set()
    # Loop through input file
    with open(results) as infile:
        if header:
            infile.next()
        for line in infile:
            # Extact and store input gene
            lineData = line.strip().split('\t')
            gene = lineData[geneCol]
            if gene in allGenes:
                raise ValueError('Genes duplicated in results file')
            # Extract and store statistic
            stat = lineData[statCol]
            if stat == 'NA':
                continue
            stat = float(lineData[statCol])
            # Store genes and signifcant genes
            allGenes.add(gene)
            if stat <= statMax:
                sigGenes.add(gene)
    # Return data
    return((allGenes, sigGenes))

def extract_gene_results_posneg(
        results, geneCol, log2Col, statCol, statMax, header
    ):
    ''' Parses text files to extract genes for ontology analysis.
    
    Args:
         results - Full path to text file.
         geneCol (int) - Column containing gene names.
         log2Col (int) - Column containing log2 values/
         statCol (int) - Column containing significance statistic.
         statMax (float) - Threshold of significance for statistic.
         header (bool) - Whether text file contains a header.
    
    Returns:
        allGenes (set)- The set of all genes in file with an
            associated statistic and log fold change.
        posGenes (set)- The set of significane genes in file
            with a positive log fold change.
        negGenes (set)- The set of significane genes in file
            with a negative log fold change.
    
    '''
    # Check arguments
    if not isinstance(geneCol, int) or geneCol < 0:
        raise ValueError('Unacceptable value for geneCol')
    if not isinstance(log2Col, int) or log2Col < 0:
        raise ValueError('Unacceptable value for geneCol')
    if not isinstance(statCol, int) or statCol < 0:
        raise ValueError('Unacceptable value for statCol')
    if not isinstance(statMax, float) or 0 > statMax > 1:
        raise ValueError('Unacceptable value for statMax')
    if not isinstance(header, bool):
        raise ValueError('Unacceptable value for header')
    # Create variables to store data
    allGenes = set()
    posGenes = set()
    negGenes = set()
    # Loop through input file
    with open(results) as infile:
        if header:
            infile.next()
        for line in infile:
            # Extact and store input gene
            lineData = line.strip().split('\t')
            gene = lineData[geneCol]
            if gene in allGenes:
                raise ValueError('Genes duplicated in results file')
            # Extract and store statistic
            stat = lineData[statCol]
            if stat == 'NA':
                continue
            stat = float(stat)
            log2 = lineData[log2Col]
            if log2 == 'NA':
                continue
            log2 = float(log2)
            # Store genes and signifcant genes
            allGenes.add(gene)
            if stat <= statMax:
                if log2 > 0:
                    posGenes.add(gene)
                elif log2 < 0:
                    negGenes.add(gene)
                else:
                    raise ValueError('Unexpected log2 value')
    # Return data
    return((allGenes, posGenes, negGenes))

def parse_gmt(gmt):
    ''' Parses Broad gmt files for gene ontology analysis.
    
    Args:
        gmt - Full path to gmt file

    Returns:
        anno2gene - A dictionary where keys are annotations
        and values are sets of genes.
    
    '''
    # Create dictionaries to store data
    anno2gene = {}
    # Create dictionary to store annotation and loop through file
    with open(gmt) as infile:
        for line in infile:
            # Extract data from line
            lineData = line.strip().split('\t')
            annotation = lineData[0]
            geneSet = set(lineData[2:])
            # Add data to dictionary
            if annotation in anno2gene:
                if not anno2gene[annotation] == geneSet:
                    raise ValueError(
                        'Conflicting gene sets: {}'.format(annotation))
            anno2gene[annotation] = geneSet
    # Return dictionary
    return(anno2gene)

def calculate_hypergeo(
        allGenes, sigGenes, geneAnno, minGO = 5, maxGO = 500, minGene = 3,
        annoGenesOnly = True
    ):
    ''' Generates pvalues and fdr for gene ontology terms.
    
    Args:
        allGenes (set)- A set of all genes in experiment
        sigGenes (set)- A set of significant genes in experiment.
        geneAnno (dict)- A dictionary of gene ontology terms and sets of
            genes associated with the term.
        minGO (int)- Minimum number of genes associated with GO term.
        maxGO (int)- Maximum number of genes associated with GO term.
        minGene (int)- Minimum number of significant genes assocaited
            with GO term.
        annoGenesOnly (bool)- Whether to only use genes with annotation
            in pvalue calculation. False will use all genes.

    Returns:
        outDF - A pandas dataframe containing the results of the analyis.
    
    '''
    # Check arguments
    if not isinstance(minGO, int) or minGO < 1:
        raise ValueError
    if not isinstance(maxGO, int) or maxGO < minGO:
        raise ValueError
    if not isinstance(minGene, int) or minGene < 1:
        raise ValueError
    # Check gene lists
    if not isinstance(allGenes, set):
        allGenes = set(allGenes)
    if not isinstance(sigGenes, set):
        sigGenes = set(sigGenes)
    if not sigGenes.issubset(allGenes):
        raise ValueError('sigGenes must be subset of allGenes')
    # Extract background genes for pvalue calculation and count
    if annoGenesOnly:
        annoGenes = set()
        for genes in geneAnno.itervalues():
            for gene in genes:
                annoGenes.add(gene)
        background = annoGenes.intersection(allGenes)
        backgroundSig = annoGenes.intersection(sigGenes)
    else:
        background = allGenes
        backgroundSig = sigGenes
    N = len(background)
    n = len(backgroundSig)
    # Loop through gene annotation
    outList = []
    for anno, annoList in geneAnno.iteritems():
        annoSet = set(annoList)
        # Check number of genes associated with term
        K = len(background.intersection(annoSet))
        if minGO > K > maxGO:
            continue
        # Check number of significant genes associated with term
        k = len(backgroundSig.intersection(annoSet))
        if k < minGene:
            continue
        # Generate and store pvalue
        pvalue = 1 - ssh.cdf(k-1, N, K, n)
        outList.append((anno, N, n, K, k, pvalue))
    # Process and return data
    outDF = pd.DataFrame(outList, columns=['term','N','n','K','k','pval'])
    outDF = outDF.sort_values(by='pval')
    outDF['fdr'] = ssm.multipletests(outDF['pval'], method='fdr_bh')[1]
    return(outDF)

def calculate_hypergeo_posneg(
        allGenes, posGenes, negGenes, geneAnno, minGO = 5, maxGO = 500,
        minGene = 3, combined = False, annoGenesOnly = False
    ):
    ''' Generates pvalues and fdr for gene ontology terms. Considers
    signficant genes with a positive and negative fold change seperately.
    
    Args:
        allGenes (set)- Set of all genes in experiment
        posGenes (set)- Set of significant genes with positive fold change.
        negGenes (set)- Set of significant genes with negative fold change
        geneAnno (dict)- A dictionary of gene ontology terms and sets of
            genes associated with the term.
        minGO (int)- Minimum number of genes associated with GO term.
        maxGO (int)- Maximum number of genes associated with GO term.
        minGene (int)- Minimum number of significant genes assocaited
            with GO term.
        annoGenesOnly (bool)- Whether to only use genes with annotation
            in pvalue calculation. False will use all genes.
    
    Returns:
        outDF - A pandas dataframe containing the results of the analyis.
    
    '''
    # Check arguments
    if not isinstance(minGO, int) or minGO < 1:
        raise ValueError
    if not isinstance(maxGO, int) or maxGO < minGO:
        raise ValueError
    if not isinstance(minGene, int) or minGene < 1:
        raise ValueError
    # Check gene lists
    if not isinstance(allGenes, set):
        allGenes = set(allGenes)
    if not isinstance(posGenes, set):
        posGenes = set(posGenes)
    if not isinstance(negGenes, set):
        negGenes = set(negGenes)
    if not posGenes.issubset(allGenes):
        raise ValueError('posGenes must be subset of allGenes')
    if not negGenes.issubset(allGenes):
        raise ValueError('negGenes must be subset of allGenes')
    if len(posGenes.intersection(negGenes)) > 0:
        raise ValueError('overlap between posGenes and negGenes')
    # Extract background genes for pvalue calculation and count
    if annoGenesOnly:
        annoGenes = set()
        for genes in geneAnno.itervalues():
            for gene in genes:
                annoGenes.add(gene)
        background = annoGenes.intersection(allGenes)
        sigDict = {'pos':annoGenes.intersection(posGenes),
            'neg':annoGenes.intersection(negGenes)}
    else:
        background = allGenes
        sigDict = {'pos':posGenes, 'neg':negGenes}
    # Add combined set if required
    if combined:
        sigDict['com'] = sigDict['pos'].union(sigDict['neg'])
    # Loop through conditions and gene annotation
    outList = []
    N = len(background)
    for condition in sigDict.keys():
        sigGenes = sigDict[condition]
        n = len(sigGenes)
        for anno, annoList in geneAnno.iteritems():
            annoSet = set(annoList)
            # Check number of genes associated with term
            K = len(background.intersection(annoSet))
            # Check number of significant genes associated with term
            k = len(sigGenes.intersection(annoSet))
            if k < minGene or minGO > K > maxGO:
                pvalue = np.NaN
            else:
                pvalue = ssh.sf(k, N, K, n, loc=1)
            outList.append((anno, condition, N, n, K, k, pvalue))
    # Add false discovery rate
    outDF = pd.DataFrame(outList, columns=[
        'term', 'query', 'N','n','K','k','pval'])
    outDF = outDF.sort_values(by='pval')
    pvalues = outDF['pval'][~np.isnan(outDF['pval'])]
    pvalueIndices = outDF['pval'].index[~outDF['pval'].apply(np.isnan)]
    fdr =  ssm.multipletests(pvalues, method='fdr_bh')[1]
    outDF.loc[pvalueIndices, 'fdr'] = fdr
    return(outDF)

def extract_overlap_results_posneg(
        outPrefix, gmt, results, geneCol, log2Col, statCol, statMax = 0.05
    ):
    # Parse gmt file
    gmtData = parse_gmt(gmt)
    # Create output dictionary
    outputDict = {}
    for term in gmtData.keys():
        outputDict[(term, 'pos')] = []
        outputDict[(term, 'neg')] = []
    # Open results file and extract header
    with open(results) as inFile:
        # Add header to output dictionary
        header = inFile.next().strip()
        for term in outputDict:
            outputDict[term].append(header)
        # Extract line data
        for line in inFile:
            lineData = line.strip().split('\t')
            # Extract and check stat data
            stat = lineData[statCol]
            if stat == 'NA':
                continue
            stat = float(stat)
            # Loop through gmt terms and find matches
            log2 = lineData[log2Col]
            if log2 == 'NA':
                continue
            log2 = float(log2)
            # Store genes and signifcant genes
            gene = lineData[geneCol]
            for term, termGenes in gmtData.items():
                if gene in termGenes:
                    if log2 > 0:
                        outputDict[(term, 'pos')].append(line.strip())
                    elif log2 < 0:
                        outputDict[(term, 'neg')].append(line.strip())
    # Loop through data and create output files
    for term, change in outputDict.keys():
        outLines = '\n'.join(outputDict[(term, change)])
        term = re.sub('\s', '_', term)
        outFile = '{}.{}.{}.results'.format(outPrefix, term, change)
        with open(outFile, 'w') as out:
            out.write(outLines)

extract_overlap_results_posneg(
    outPrefix = '/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/geneOntology/customResults/overlapData/A_WTvsMU',
    gmt = '/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/geneOntology/customGmt/myrto.guilherme.ensembl.gmt',
    results = '/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/rnaSeqPool/myrto/DESeq/A_WTvsMU.results',
    geneCol = 0,
    log2Col = 3,
    statCol = 7
)

def extract_ensembl_names(gtf):
    # Create regular expressions
    eRE = re.compile('gene_id\s+"(.*?)";')
    nRE = re.compile('gene_name\s+"(.*?)";')
    # Sequentially process lines in file and add to dictionary
    nameDict = collections.defaultdict(set)
    with open(gtf) as filein:
        for line in filein:
            # Skip lines starting with #
            if line.startswith('#'):
                continue
            # Extract ensembl and gene names from dictionary
            data = line.strip().split('\t')[8]
            ensembl = re.search(eRE, data).group(1)
            name = re.search(nRE, data).group(1)
            nameDict[ensembl].add(name)
    # Convert sets to lists
    nameDict = {k:list(v) for k,v in nameDict.iteritems()}
    # Create log and return data
    logDict = collections.defaultdict(int)
    for values in nameDict.itervalues():
        logDict[len(values)] += 1
    return(nameDict, logDict)

def span_translations(tran1, tran2, log=True):
    # Check argument and create output variables
    if not isinstance(log, bool):
        raise TypeError('Log argument must be boolean')
    tranDict = {}
    # Loop through values in 1st dictionary
    for key, valueList1 in tran1.iteritems():
        # Create entry for key in output dictionary
        if key not in tranDict:
            tranDict[key] = set()
        # Loop through values in 1st dictionary
        for v1 in valueList1:
            # Extract and stroe values for entries in 2nd dictionary
            if v1 in tran2:
                valueList2 = tran2[v1]
                for v2 in valueList2:
                    tranDict[key].add(v2)
    # Convert sets to lists
    tranDict = {k:list(v) for k,v in tranDict.iteritems()}
    # Create log and return data
    logDict = collections.defaultdict(int)
    for values in tranDict.itervalues():
        logDict[len(values)] += 1
    return(tranDict, logDict)






#parse_gmt('/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/genomeData/geneOntology/c2.all.v5.1.symbols.gmt')
#
#
#ensemblTran, ensemblLog = create_tran_dictionary('/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/genomeData/mouse2Human/ensembl_mouse_human.txt')
#nameTran, nameLog = extract_ensembl_names('/farm/scratch/rs-bio-lif/rabino01/myrtoDenaxa/genomeData/mouse2Human/Homo_sapiens.GRCh38.84.gtf')
#finalTran, finalLog = span_translations(ensemblTran, nameTran)
#for inGene, outGene in finalTran.iteritems():
#    if len(outGene) == 0:
#        print('{}\t'.format(inGene))
#    else:
#        for gene in outGene:
#            print('{}\t{}'.format(inGene,gene))
