"""HiC_BinCounts.py

Usage:
    
    HiC_BinCounts.py bed <bedfile> <mincount> <outdir> <inputfiles>...
        [--label=<label>] [--threads=<threads>]
    
    HiC_BinCounts.py nobed <chrfile> <binsize> <mincount> <outdir>
        <inputfiles>... [--label=<label>] [--threads=<threads>] [--equal]
    
    HiC_BinCounts.py (-h | --help)
    
Options:
    
    --label=<label>      Label to add output file names
    --threads=<threads>  Number of threads [default: 1]
    --equal              Bins should be equally sized
    --help               Output this message
    
"""
# Import required modules
import os
import re
import numpy as np
from ngs_python.structure import interactionMatrix, analyseInteraction
from general_python import docopt, toolbox
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Check numerical arguments
args['<mincount>'] = int(args['<mincount>'])
args['--threads'] = int(args['--threads'])
if args['nobed']:
    args['<binsize>'] = int(args['<binsize>'])
# Check input files and output directory
if args['bed']:
    toolbox.checkArg(args['<bedfile>'], 'file')
else:
    toolbox.checkArg(args['<chrfile>'], 'file')
for f in args['<inputfiles>']:
    toolbox.checkArg(f, 'file')
toolbox.checkArg(args['<outdir>'], 'dir')
# Modify label argument
if args['--label']:
    toolbox.checkArg(args['--label'], 'str')
    args['--label'] = '_' + args['--label']
# Extract and print sample names
sampleNames = [re.search('([^/]*)\.fragLigations\.gz$',f).group(1) for f in args['<inputfiles>']]
logData = 'Samples:\n  %s\n' %(
    '\n  '.join(sampleNames) 
)
# Extract and print parameters to create bins
if args['bed']:
    binData = args['<bedfile>']
    logData += '\nParameters:\n\t%s\n\t%s\n' %(
        'bed file provided',
        'minimum bin count: %s' %(args['<mincount>'])
    )
else:
    binData = (args['<chrfile>'], args['<binsize>'], args['--equal'])
    logData += '\nParameters:\n  %s\n  %s\n  %s\n' %(
        'max bin size: %s' %(args['<binsize>']),
        'bin size equal: %s' %(args['--equal']),
        'minimum bin count: %s' %(args['<mincount>'])
    )
# Create bin object and save bed
genomeBins = interactionMatrix.genomeBin(binData)
# Sequentially process input files
failedBins = np.full(genomeBins.binCount, False, dtype=bool)
prefixList = []
for f in args['<inputfiles>']:
    # Extract sample names
    sampleName = re.search('([^/]*)\.fragLigations.gz$',f).group(1)
    # Create output file prefix
    if args['--label']:
        prefix = os.path.join(args['<outdir>'], sampleName + args['--label'])
    else:
        prefix = os.path.join(args['<outdir>'], sampleName
    prefixList.append(prefix)
    # Create interaction matrix and save to file
    countMatrix, logArray = genomeBins.generateMatrix(
        f, args['--threads'])
    np.savetxt(prefix + '.countMatrix.gz', countMatrix, '%s', '\t',
        header = '\t'.join(genomeBins.binNames), comments = '')
    # Print interaction data
    logData += '\n%s:\n  Interaction Data:\n%s\n%s\n%s\n%s\n' %(
        sampleName,
        '    total: %s' %(logArray[0]),
        '    accepted: %s' %(logArray[3]),
        '    no chromosome: %s' %(logArray[1]),
        '    no bin: %s' %(logArray[2])
    )
    # Find and process bins below minimum bin count
    colSums = countMatrix.sum(axis=0)
    binsBelowMin = colSums < args['<mincount>']
    failedBins = np.logical_or(failedBins, binsBelowMin)
    # Store bin data
    logData += '  Bin Data:\n%s\n%s\n%s\n%s\n' %(
        '    bin number: %s' %(genomeBins.binCount),
        '    bins below minimum: %s' %(sum(binsBelowMin)),
        '    mean bincount: %.0f' %(np.mean(colSums)),
        '    median bin count: %.0f' %(np.median(colSums))
    )
# Extract combined bin data
logData += '\nCombined Samples:\n  Bin Data:\n%s\n%s\n' %(
    '    bin number: %s' %(genomeBins.binCount),
    '    bins below minimum: %s' %(sum(failedBins))
)
# Save log files
for prefix in prefixList:
    with open(prefix + '.matrixLog', 'w') as outFile:
        outFile.write(logData)
# Normalise matrices
if sum(failedBins) < genomeBins.binCount:
    # Loop through count matrices
    for prefix in prefixList:
        # Normalise matrix
        normMatrix, biasData = interactionMatrix.normaliseMatrix(
            prefix + '.countMatrix.gz', failedBins)
        # Save normalised interactions
        np.savetxt(prefix + '.normMatrix.gz', normMatrix, '%.6f', '\t',
            header = '\t'.join(genomeBins.binNames), comments = '')
        # Extract bin level data
        maskMatrix = analyseInteraction.maskMatrix(prefix + '.normMatrix.gz')
        maskMatrix.binDF['bias'] = biasData
        maskMatrix.binDirection()
        maskMatrix.binDistance()
        maskMatrix.binDF.to_csv(prefix + '.binData', '\t', '')
        # Extract global data
        distance = maskMatrix.combinedDistance(0.1)
        np.savetxt(prefix + '.dist.gz', distance, '%s', '\t')
