"""resultsGOAnalysis.py

Usage:
    
    resultsGOAnalysis.py <results> <geneCol> <statCol> <statMax> <gmt>
        <outFile> [--minGO=<minGO>] [--maxGO=<maxGO>] [--minGene=<minGene>]
        [--log2Col=<log2Col>] [--includeCombined] [--onlyAnno] [--noHeader]
    
Options:
    
    --maxGO=<maxGO>      Maximum genes in GO set [default: 1000].
    --minGO=<minGO>      Minimum genes in GO set [default: 5].
    --minGene=<minGene>  Minimum significant genes in GO set [default: 3].
    --log2Col=<log2Col>  Column for log2 fold change data. Supplying this
        results in positive and negative fold change genes being considered
        seperately.
    --includeCombined    Include combined geneset alongside positive and
        negative genes sets. Only effectice with --log2Col argument.
    --onlyAnno           Only consider genes with annotation.
    --noHeader           Results file has no header.
    
"""
# Import required modules
from ngs_python.gtf import gene_conversion
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<geneCol>'] = int(args['<geneCol>'])
args['<statCol>'] = int(args['<statCol>'])
args['<statMax>'] = float(args['<statMax>'])
args['--minGO'] = int(args['--minGO'])
args['--maxGO'] = int(args['--maxGO'])
if args['--log2Col'] is not None:
    args['--log2Col'] = int(args['--log2Col'])
# Parse gmt
geneAnno = gene_conversion.parse_gmt(args['<gmt>'])
if isinstance(args['--log2Col'], int):
    # Extract gene list
    allGenes, posGenes, negGenes = gene_conversion.extract_gene_results_posneg(
        results=args['<results>'], geneCol=args['<geneCol>'],
        log2Col=args['--log2Col'], statCol=args['<statCol>'],
        statMax=args['<statMax>'], header=not(args['--noHeader'])
    )
    # Perform annotation
    goData = gene_conversion.calculate_hypergeo_posneg(
        allGenes=allGenes, posGenes=posGenes, negGenes=negGenes,
        geneAnno=geneAnno, annoGenesOnly=args['--onlyAnno'],
        combined=args['--includeCombined']
    )
else:
    # Extract gene list
    allGenes, sigGenes = gene_conversion.extract_gene_results(
        results=args['<results>'], geneCol=args['<geneCol>'],
        statCol=args['<statCol>'], statMax=args['<statMax>'],
        header=not(args['--noHeader'])
    )
    # Perform annotation
    goData = gene_conversion.calculate_hypergeo(
        allGenes=allGenes, sigGenes=sigGenes, geneAnno=geneAnno,
        annoGenesOnly=args['--onlyAnno']
    )
print(goData)
goData.to_csv(args['<outFile>'], sep='\t', index=False)
