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
from ngs_analysis.structure import interactionMatrix
from general_functions import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Check numerical arguments
args['<mincount>'] = int(args['<mincount>'])
args['--threads'] = int(args['--threads'])
if args['nobed']:
    args['<binsize>'] = int(args['<binsize>'])
# Check input files and output directory
if args['bed']:
    if not os.path.isfile(args['<bedfile>']):
        raise IOError('Bed file %s not found' %(args['<bedfile>']))
else:
    if not os.path.isfile(args['<chrfile>']):
        raise IOError('Chromosome file %s not found' %(args['<chrFile>']))
for f in args['<inputfiles>']:
    if not os.path.isfile(f):
        raise IOError('Input file %s not foundi' %(f))
if not os.path.isdir(args['<outdir>']):
    raise IOError('Output directory %s not found' %(args['<outdir>']))
# Extract and print parameters to create bins
if args['bed']:
    binData = args['<bedfile>']
    print 'Parameters:\n\t%s\n\t%s' %(
        'bed file provided',
        'minimum bin count: %s' %(args['<mincount>'])
    )
else:
    binData = (args['<chrfile>'], args['<binsize>'], args['--equal'])
    print 'Parameters:\n  %s\n  %s\n  %s' %(
        'max bin size: %s' %(args['<binsize>']),
        'bin size equal: %s' %(args['--equal']),
        'minimum bin count: %s' %(args['<mincount>'])
    )
# Create bin object and save bed
genomeBins = interactionMatrix.genomeBin(binData)
# Sequentially process input files
failedBins = np.full(genomeBins.binCount, False, dtype=bool)
matrixFileList = []
for f in args['<inputfiles>']:
    # Extract sample names
    sampleName = re.search('([^/]*)\.fragLigations.gz$',f).group(1)
    # Create output file prefix
    if args['--label']:
        outPrefix = args['<outdir>'] + sampleName + '_' + args['--label']
    else:
        outPrefix = args['<outdir>'] + sampleName
    # Create output file names
    bedFile = outPrefix + '.bed'
    matrixFile = outPrefix + '.countMatrix.gz'
    biasFile = outPrefix + '.bias'
    normMatrixFile = outPrefix + '.normMatrix.gz'
    # Store file names required for normalisation
    matrixFileList.append((matrixFile, biasFile, normMatrixFile))
    # Save bed file
    genomeBins.writeBed(bedFile)
    # Create interaction matrix and save to file
    countMatrix, logArray = interactionMatrix.generateMatrix(
        f, genomeBins, args['--threads'])
    np.savetxt(matrixFile, countMatrix, '%s', '\t',
        header = '\t'.join(genomeBins.binNames), comments = '')
    # Print Interaction Data
    print '\n%s:\n  Interaction Data:\n%s\n%s\n%s\n%s' %(
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
    # Print bin data
    print '  Bin Data:\n%s\n%s\n%s\n%s' %(
        '    bin number: %s' %(genomeBins.binCount),
        '    bins below minimum: %s' %(sum(binsBelowMin)),
        '    mean bincount: %.0f' %(np.mean(colSums)),
        '    median bin count: %.0f' %(np.median(colSums))
    )
# Print combined bin data
print '\nCombined Samples:\n  Bin Data:\n%s\n%s' %(
    '    bin number: %s' %(genomeBins.binCount),
    '    bins below minimum: %s' %(sum(failedBins))
)
# Normalise matrices
if sum(failedBins) < genomeBins.binCount:
    # Loop through count matrices
    for files in matrixFileList:
        # Normalise matrix
        normMatrix, biasData = interactionMatrix.normaliseMatrix(
            files[0], failedBins)
        # Save interactions
        np.savetxt(files[1], np.array([genomeBins.binNames,biasData]).T,
            '%s', '\t')
        np.savetxt(files[2], normMatrix, '%.6f', '\t',
            header = '\t'.join(genomeBins.binNames), comments = '')
