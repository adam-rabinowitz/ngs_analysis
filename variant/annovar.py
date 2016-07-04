from ngs_python.variant import varscan
import pandas as pd
import os
import subprocess
import numpy as np

def geneAnno(
        inFile, outPrefix, path, buildver, database
    ):
    '''
    Function perform gene annovar gene annotation on a valid
    annovar input file. Function takes X arguments
    
    '''
    # Build command
    command = [path, '-geneanno', '-buildver', buildver, '-outfile',
        outPrefix, inFile, database]
    # Join and return command
    command = ' '.join(command)
    return(command)

def geneAnno2DF(
        variantList, path, buildver, database, tempprefix
    ):
    # Create temporary file names
    annoIn = tempprefix + '.av'
    varFunc = tempprefix + '.variant_function'
    exonVarFunc = tempprefix + '.exonic_variant_function'
    annoLog = tempprefix + '.log'
    invalidIn = tempprefix + ''
    # Create gene annotation input
    with open(annoIn, 'w') as outFile:
        for chrom, start, reference, variant in variantList:
            # Create variant name
            variantName = '%s:%s:%s:%s' %(
                chrom, start, reference, variant
            )
            # Process deletions
            if "-" in variant:
                # Set end and reference values
                end = start + (len(variant) - 1)
                reference = '0'
                variant = '-'
            else:
                end = start
            # Write variant data to output file
            outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(chrom, start, end,
                reference, variant, variantName))
    # Perform annotation
    annoCommand = geneAnno(inFile = annoIn, outPrefix = tempprefix,
        path = path, buildver = buildver, database = database)
    subprocess.check_call(annoCommand, shell = True)
    # Read in annovar output files
    vf = pd.read_csv(varFunc, sep = "\t", index_col = -1, header = None)
    vf.columns = ['class', 'genes', 'chrom', 'start', 'end', 'ref', 'var']
    # Create output dataframe from annovar output
    try:
        evf = pd.read_csv(exonVarFunc, sep = "\t", index_col = -1,
            header = None)
    except ValueError:
        evf = None
    if evf is None:
        outDF = vf.loc[:,['class','genes']]
        outDF['affect'] = np.nan
        outDF['change'] = np.nan
    else:
        evf.columns = ['line', 'affect', 'change', 'chrom', 'start', 'end',
            'ref', 'var']
        outDF = pd.concat(
            [vf[['class','genes']], evf[['affect','change']]],
            axis = 1
        )
    # Delete temporary files
    for f in [annoIn, varFunc, exonVarFunc, annoLog]:
        if os.path.isfile(f):
            os.remove(f)
    # Return annovar output
    return(outDF)
