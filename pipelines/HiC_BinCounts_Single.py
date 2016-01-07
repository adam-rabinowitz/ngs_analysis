# Import required modules
from argparse import ArgumentParser
import numpy as np
from ngs_analysis.structure import interactionMatrix 
# Create argument parser
parser = ArgumentParser()
# List required arguments
parser.add_argument('inFile', help = 'File listing input fragend ligations',
    type = str)
parser.add_argument("outPrefix", help = "Prefix of output file name",
    type = str)
parser.add_argument('-c', '--chrSize', help = 'Tab delimited file listing '+\
    ' chromosome names and sizes', default = '', type = str)
parser.add_argument('-s', '--binSize', help = 'Maximum size of genomic bins',
    default = 0, type = int)
parser.add_argument('-e', '--equal', help = 'Equal bin size', default = True,
    type = bool)
parser.add_argument('-m', '--minCount', help = 'Minimum count for a bin to be '+\
    'normalised', default = 0, type = int)
parser.add_argument('-b', '--bedFile', help = 'bed file containing predefined '+\
    'bins', default = '', type = str)
parser.add_argument('-t', '--threads', help = 'number of threads to use',
    default = 1, type = int)
# Extract arguments
args = parser.parse_args()
# Check input files and print parameters
if args.bedFile:
    binData = args.bedFile
    print 'Parameters:\n\t%s\n\t%s' %(
        'bed file provided',
        'minimum bin count: %s' %(args.minCount)
    )
elif args.chrSize and args.binSize:
    binData = (args.chrSize, args.binSize, args.equal)
    print 'Parameters:\n\t%s\n\t%s\n\t%s' %(
        'max bin size: %s' %(args.binSize),
        'bin size equal: %s' %(args.equal),
        'minimum bin count: %s' %(args.minCount)
    )
else:
    raise IOError('bed file or chromosome file and bin size must be provided')
# Create output files
args.bedFile = args.outPrefix + '.bed'
args.biasFile = args.outPrefix + '.bias'
args.countMatrix = args.outPrefix + '.countMatrix'
args.normMatrix = args.outPrefix + '.normMatrix'
# Create bin object and save bed
genomeBins = interactionMatrix.genomeBin(binData)
genomeBins.writeBed(args.bedFile)
# Count interactions
countMatrix, logArray = interactionMatrix.generateMatrix(
    args.inFile, genomeBins, args.threads
)
# Save interactions
np.savetxt(args.countMatrix, countMatrix, '%s', '\t',
    header = '\t'.join(genomeBins.binNames), comments = '')
# Print Interaction Data
print '\nInteraction Data\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'total: %s' %(logArray[0]),
    'accepted: %s' %(logArray[3]),
    'no chromosome: %s' %(logArray[1]),
    'no bin: %s' %(logArray[2])
)
# Count number of bins above minimum
colSums = countMatrix.sum(axis=0)
binsAboveMin = sum(colSums >= args.minCount)
# Print bin data
print '\nBin Data\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'total: %s' %(len(colSums)),
    'bins above minimum: %s' %(binsAboveMin),
    'mean count: %.0f' %(np.mean(colSums)),
    'median count: %.0f' %(np.median(colSums))
)
if binsAboveMin > 0:
    # Normalise matrix
    normMatrix, biasData = interactionMatrix.normaliseMatrix(
        countMatrix, args.minCount)
    # Save interactions
    np.savetxt(args.normMatrix, normMatrix, '%.6f', '\t',
        header = '\t'.join(genomeBins.binNames), comments = '')
    np.savetxt(args.biasFile, np.array([genomeBins.binNames,biasData]).T,
        '%s', '\t')
    
