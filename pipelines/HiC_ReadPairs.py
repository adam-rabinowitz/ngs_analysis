###############################################################################
## Initialise analysis
###############################################################################
# Import required modules
import subprocess
import sys
import numpy
import argparse
import itertools
import os
# Import custom modules
from ngs_analysis.fastq import fastqFind, fastqMerge, fastqAlign
from ngs_analysis.structure import alignedPairs, fragendPairs
# Create parser
parser = argparse.ArgumentParser()
# Define positional arguments
parser.add_argument('sampleData', help = 'Comma seperated list of FASTQ '+\
    'prefix and sample name', type = str)
parser.add_argument('outDir', help = 'Output directory', type = str)
parser.add_argument('cutSite', help = 'Restriction enzyme recognition site',
    type = str)
parser.add_argument('fastqDir', help = 'Comma seperated list of directories '+\
    'containing FASTQ files', type = str)
parser.add_argument("bwaFasta", help = 'BWA indexed genome FASTA file',
    type = str)
# Define optional arguments
parser.add_argument('-m', '--minMapQ', help = 'Minimum read mapping quality '+\
    '(Minimum = 2)', type = int, default = 10)
parser.add_argument('-l', '--minLength', help = "Minimum length of trimmed '+\
    'reads (Minimum = 18)", type = int, default = 20)
parser.add_argument('-d', '--maxDistance', help = 'Max distance (bp) from '+\
    'read to restriction site', type = int, default = 1000)
parser.add_argument('-c', '--maxConcordant', help = 'Max fragment size of '+\
    'concordant reads', type = int, default = 2000)
parser.add_argument('-r', '--removeConcordant', help = 'Remove concordant '+\
    'pairs ', type = bool, default = True)
# Check arguments
args = parser.parse_args()
if args.maxDistance < 0:
    exit("Script terminated: -d option must be positive")
if args.minMapQ < 2:
    exit("Script terminated: -m option must be 2 or higher")
if args.minLength < 18:
    exit("Script terminated: -l option must be 18 or higher")
# Process areguments
args.cutSite = args.cutSite.upper()
args.fastqPrefix, args.sampleName = args.sampleData.split(',')
# Print parameters
print 'Parameters:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'minimum read length: ' + str(args.minLength)
    'mimimum mapping quality: ' + str(args.minMapQ)
    'maximum fragend distance: ' + str(args.maxDistance)
    'maximum distance of concordant pairs: ' str(args.maxConcordant)
    'remove concordant pairs: ' str(args.removeConcordant)
)
# Create output file names
args.logFile = args.outDir + args.sampleName + '.log'
args.outFastq = args.outDir + args.sampleName + '_trimmed.fastq.gz'
args.outSam = args.outDir + args.sampleName + '.sam'
args.nameSortBam = args.outDir + args.sampleName + "_nSort.bam"
args.outPairs = args.outDir + args.sampleName + ".readPairs.gz"
args.outFrags = args.outDir + args.sampleName + ".fragLigations.gz"

################################################################################
## Process FASTQ files and perform alignment
################################################################################
# Extract fastq file names
args.read1, args.read2 = fastqFind.findIlluminaFastq(prefix = args.fastqPrefix,
    dirList = args.fastqDir.split(','), pair = True)
# Trim and merge fastq files
trimMetrics = fastqMerge.mergeLabelTrimPair(
    fastqIn1 = args.read1,
    fastqIn2 = args.read2,
    fastqOut = args.outFastq,
    args.trimSeq = args.cutSite,
    label1 = ':1',
    label2 = ':2'
)
# Print trim metrics
print 'Trim Metrics:\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'total: ' + str(trimMetrics['total']),
    'too short: ' + str(trimMetrics['short']),
    'read1 trim: ' + str(trimMetrics['trim1']),
    'read2 trim: ' + str(trimMetrics['trim1'])
)
# Generate align command
alignCommand = fastqAlign.bwaMemAlign(
    index = args.bwaFasta,
    outSam = args.outSam,
    read1 = args.outFastq,
    path = 'bwa',¬
    threads = 4,
    markSecondary = True,
    check = True
)
# Generate sort command
sortCommand = samtools.sort(
    inFile = args.outSam,
    outFile = args.nameSortBam,
    name = True,
    threads = 4,¬
    memory = '2',
    delete = True,
    path = 'samtools'
)
# Merge commands and run
alignSort = '%s && %s' (alignCommand, sortCommand)
subprocess.check_output(alignSort, shell = True, stderr=subprocess.STDOUT)

###############################################################################
## Extract aligned pairs
###############################################################################
# Create pair object and extract pairs
pairObject = alignedPairs.ReadPairs(args.nameSortBam)
alignData = pairObect.extract(
    outPairs = args.outPairs,
    minMapQ = 20,
    maxSize = 2000,
    rmConcord = True,
    rmDup = True
)
# Print alignment metrics
print 'Alignment Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'total: ' + str(alignData['total']),
    'unmapped: ' + str(alignData['unmapped']),
    'poorly mapped: ' + str(alignData['poormap']),
    'secondary alignment: ' + str(alignData['secondary']),
    'singletons: ' + str(alignData['singletons']),
    'multiple alignments: ' + str(alignData['multiple']),
    'paired: %s' + str(alignData['pairs'])
)
# Print pair statistics
print 'Pair Data:\n\t%s\n\t%s\n\t%s\n\t''\ttotal: %s' %(alignData['pairs'] / 2)
print '\tunique: %s' %(alignData['unique'])
print '\tduplication rate: %s' %(alignData['duprate'])
print '\tdiscordant: %s' %(alignData['discord'])
print '\tconcordant rate %s' %(alignData['concordrate'])

##############################################################################
## Extract fragend pairs
##############################################################################
fragendObject = fragendPairs.Fragend(args.bwaFasta, args.cutSite)

# Remove fragend and ligation distance data
fragendDist = fragendData.pop('fragend distance')
ligationDist = fragendData.pop('ligation distance')
# Print fragend metrics
print '\nFragend Data:'
for key in fragendData:
    print '\t%s: %s' %(key, fragendData[key])
# Calclate and print fragend metrics if accepted fragends are present
fragendCount = float(len(fragendDist))
if fragendCount > 0:
    fragendHist = numpy.histogram(
        fragendDist,
        bins = [0,100,200,300,400,500,600,700,800,900,1000,float('inf')]
    )
    fragendCumFreq = numpy.cumsum(fragendHist[0]) / float(fragendCount)
    # Print fragend distance data
    print '\nFragend Distance:'
    for bin, freq in itertools.izip(fragendHist[1][1:], fragendCumFreq):
        print '\t<%.0f: %.3f' %(
            bin,
            freq
        )
    # Process ligation data
    print '\nIntrachromosomal Ligation Data:'
    print ('\tmean distance: %.0f\n\t25th percentile: %.0f\n\t50th percentile:' 
        '%.0f\n\t75th percentile: %.0f') %(
        numpy.mean(ligationDist),
        numpy.percentile(ligationDist,25),
        numpy.percentile(ligationDist,50),
        numpy.percentile(ligationDist,75)
    )
