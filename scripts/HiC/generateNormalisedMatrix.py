"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <mincount> <infiles>...
        [--scope=<scope>] [--threads=<threads>]
    
    generateNormalisedMatrix.py (-h | --help)
    
Options:
    
    --threads=<threads>  Number of threads [default: 1]
    --scope=<scope>      Scope of normalisation. Either genome ('gen'),
                         chromosome ('chr') or path to a tab delimited
                         text file outlining chromosome, start, end and
                         name of regions.
    --help               Output this message
    
"""
# Import required modules
import os
import re
import gzip
import numpy as np
import pandas as pd
from ngs_python.structure import interactionMatrix, analyseInteraction
from general_python import docopt, toolbox
# Extract arguments
args = docopt.docopt(__doc__, version = '1.0')
# Check bin names are identical
for count, infile in enumerate(args['<infiles>']):
    # Open input files
    with gzip.open(infile, 'r') as openFile:
        header = openFile.next().strip().split('\t')
    # Check all headers are identical
    if count:
        if not header == binNames:
            raise IOError('Input files must have identical headers')
    else:
        binNames = header
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
# Extract sample names
sampleNames = [os.path.basename(x)[:-15] for x in args['<infiles>']]
# Find bins to exclude due to low counts
excludeBins = np.zeros(len(binNames), dtype = 'bool')
for infile in args['<infiles>']:
    # Read in file
    counts = np.loadtxt(infile, dtype = np.uint32, delimiter = '\t',
        skiprows = 1)
    # Extract sums and check for symettry
    rowSums = counts.sum(axis=0)
    colSums = counts.sum(axis=1)
    if not np.all(rowSums == colSums):
        raise IOError('Count files are not symetrical')
    # Extract low bins and store
    lowBins = rowSums < np.uint32(args['<mincount>'])
    excludeBins = np.logical_or(excludeBins, lowBins)
# Find indices of regions for normalisation
regionDict = {}
# Normalise 
with open(args['--scope'], 'r') as infile:
    for line in infile:
        # Extract region data
        chrom, start, end, name = line.strip().split('\t')
        acceptableBins = (binChr == chrom) & (binStart >= np.uint32(start)) & (binEnd <= np.uint32(end))
        indices = np.where(acceptableBins)[0]
        # Add indices to dictionary
        if name in regionDict:
            regionDict[name][indices] = np.concatenate(regionDict[name]['indices'], indices)
        else:
            regionDict[name] = {'indices' : indices, 'name' : name}
print regionDict
