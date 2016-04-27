"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <regions> <mincount> <outdir> <infiles>...
        [--path=<paths>] [--threads=<threads>] [--memory=<memory>]
        [--iter=<iter>]
    
    generateNormalisedMatrix.py (-h | --help)
    
Options:
    
    --path=<path>        Path to ic_mes executable [default: ic_mes]
    --threads=<threads>  Number of threads [default: 1]
    --memory=<memory>    Memory (MB) per thread [default: 2000]
    --iter=<iter>        Iterations of ic_mes [default: 1000]
    --help               Output this message
    
"""
# Import required modules
import os
import re
import gzip
import numpy as np
import multiprocessing
from ngs_python.structure import interactionMatrix
from general_python import docopt, toolbox
# Extract arguments
args = docopt.docopt(__doc__, version = '1.0')
# Check numeric arguments
args['--threads'] = int(args['--threads'])
args['--memory'] = int(args['--memory'])
args['--iter'] = int(args['--iter'])
# Check bin names are identical
for count, infile in enumerate(args['<infiles>']):
    if not infile.endswith('countMatrix.gz'):
        raise IOError("Input files must end '.countMatrix.gz'")
    # Open input files
    with gzip.open(infile, 'r') as openFile:
        header = openFile.next().strip().split('\t')
    # Check all headers are identical
    if count:
        if not header == binNames:
            raise IOError('Input files must have identical headers')
    else:
        binNames = header
print args['<infiles>']
# Extract bin data
splitBin = [re.split('[:-]',x) for x in binNames]
binChr = np.array([x[0] for x in splitBin])
binStart = np.array([x[1] for x in splitBin], dtype = np.uint32)
binEnd = np.array([x[2] for x in splitBin], dtype = np.uint32)
# Check for overlapping bins
for chrom in np.unique(binChr):
    # Extract indices for chromsome and check they are identical
    indices = np.where(binChr == chrom)[0]
    indexDiff = np.ediff1d(indices)
    if not np.all(indexDiff == 1):
        raise IOError('Count files contains unsorted bins')
    # Extract chromosomal bin start and stop sites and interpolate
    chrStart = binStart[indices]
    chrEnd = binEnd[indices]
    binArray = np.empty(2 * len(indices), dtype = np.uint32)
    binArray[0::2] = chrStart
    binArray[1::2] = chrEnd
    # Check bins are ordered and non-overlapping
    for i in xrange(len(binArray) - 1):
        if binArray[i + 1] <= binArray[i]:
            raise IOError('Count files contain unsorted/overlapping bins')
# Create dictionaries to store indices and matrix information
indexDict ={}
matrixList = []
# Open region file and extract correpsonding bins
with open(args['<regions>'], 'r') as infile:
    for line in infile:
        # Extract region data
        chrom, start, end, rname = line.strip().split('\t')
        acceptableBins = ((binChr == chrom) & (binStart >= np.uint32(start)) &
            (binEnd <= np.uint32(end)))
        indices = np.where(acceptableBins)[0]
        # Add indices to index dictionary
        if rname in indexDict:
            indexDict[rname] = np.concatenate(indexDict[rname], indices)
        else:
            indexDict[rname] = indices
# Check index dictionary and create entries in exclude dictionary
for rname in indexDict:
    # Extract and sort indices
    indices = indexDict[rname]
    indices.sort()
    # Check for absent or duplicate indices
    if len(indices) == 0:
        raise IOError('Region %s has no corresponding bins' %(rname))
    if len(set(indices)) != len(indices):
        raise IOError('Region %s has overlapping segments' %(rname))
    # Update dictionary
    indexDict[rname] = indices
# Read in input files, save divided regions and extract counts
for rname, indices in indexDict.items():
    # Extract bin names and create list to store files
    rbinNames = np.array(binNames)[indices]
    rfileList = []
    rbinExclude = np.zeros(len(indices), dtype = 'bool')
    # Loop through files
    for infile in args['<infiles>']:
        # Extract name, create prefix
        sampleName = os.path.basename(infile)[:-15]
        outPrefix = os.path.join(args['<outdir>'],
            sampleName + '.' + rname + '.' + str(args['<mincount>']))
        # Read in file
        counts = np.loadtxt(infile, dtype = np.uint32, delimiter = '\t',
            skiprows = 1)
        # Subdivide counts, save file and store name
        regionCounts = counts[indices,:][:,indices]
        countFileName = outPrefix + '.countMatrix.gz'
        np.savetxt(countFileName, regionCounts, fmt = '%.0f', delimiter = '\t',
            header = '\t'.join(rbinNames), comments = '')
        rfileList.append(countFileName)
        # Extract sums and check for symettry
        rowSums = regionCounts.sum(axis=0)
        colSums = regionCounts.sum(axis=1)
        if not np.all(rowSums == colSums):
            raise IOError('Count files are not symetrical')
        # Extract low bins and store
        lowBins = rowSums < np.uint32(args['<mincount>'])
        rbinExclude = np.logical_or(rbinExclude, lowBins)
    # Update matrix list with commands for normalisation
    for rfile in rfileList:
        matrixList.append((rfile, rbinExclude, args['--path'], args['--memory'],
            args['--iter']))
# Generate normalised matrices using multiprocessing pool
pool = multiprocessing.Pool(args['--threads'])
pool.map(interactionMatrix.normaliseMatrixProcess, matrixList)
