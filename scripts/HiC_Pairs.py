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
from ngs_python.fastq import fastqFind, fastqMerge, fastqAlign
from ngs_python.bam import samtools
from ngs_python.structure import alignedPair, fragendPair
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
parser.add_argument('-q', '--minMapQ', help = 'Minimum read mapping quality '+\
    '(Minimum = 2)', type = int, default = 10)
parser.add_argument('-l', '--minLength', help = "Minimum length of trimmed '+\
    'reads (Minimum = 18)", type = int, default = 20)
parser.add_argument('-s', '--maxDistance', help = 'Max distance of read '+\
    'from restriction site', type = int, default = 1000)
parser.add_argument('-m', '--maxSize', help = 'Max fragment size of '+\
    'concordant reads', type = int, default = 2000)
parser.add_argument('-c', '--rmConcordant', help = 'Remove concordant pairs',
    type = bool, default = True)
parser.add_argument('-d', '--rmDuplicates', help = 'Remove duplicate pairs',
    type = bool, default = True)
parser.add_argument('-b', '--bwa', help = 'Path to BWA', type = str,
    default = 'bwa')
parser.add_argument('-t', '--threads', help = 'Number of threads to use',
    type = int, default = 4)
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
print 'Parameters:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'minimum read length: %s' %(args.minLength),
    'mimimum mapping quality: %s' %(args.minMapQ),
    'maximum fragend distance: %s' %(args.maxDistance),
    'maximum distance of concordant pairs: %s' %(args.maxSize),
    'remove duplicate pairs: %s' %(args.rmDuplicates),
    'remove concordant pairs: %s' %(args.rmConcordant)
)
# Create output file names
args.logFile = args.outDir + args.sampleName + '.log'
args.outFastq = args.outDir + args.sampleName + '_trimmed.fastq.gz'
args.outSam = args.outDir + args.sampleName + '.sam'
args.nameSortBam = args.outDir + args.sampleName + "_nSort.bam"
args.outPairs = args.outDir + args.sampleName + ".readPairs.gz"
args.outFrags = args.outDir + args.sampleName + ".fragLigations.gz"

###############################################################################
## Process FASTQ files and perform alignment
###############################################################################
# Extract fastq file names
args.read1, args.read2 = fastqFind.findIlluminaFastq(prefix = args.fastqPrefix,
    dirList = args.fastqDir.split(','), pair = True)
# Trim and merge fastq files
trimMetrics = fastqMerge.mergeLabelTrimPair(
    fastqIn1 = args.read1,
    fastqIn2 = args.read2,
    fastqOut = args.outFastq,
    trimSeq = args.cutSite,
    label1 = ':1',
    label2 = ':2',
    minLength = args.minLength
)
# Print trim metrics
print '\nTrim Metrics:\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'total: ' + str(trimMetrics['total']),
    'too short: ' + str(trimMetrics['short']),
    'read1 trim: ' + str(trimMetrics['trim1']),
    'read2 trim: ' + str(trimMetrics['trim2'])
)
# Generate align command
alignCommand = fastqAlign.bwaMemAlign(
    index = args.bwaFasta,
    outSam = args.outSam,
    read1 = args.outFastq,
    path = args.bwa,
    threads = str(args.threads),
    markSecondary = True,
    check = True
)
# Generate sort command
sortCommand = samtools.sort(
    inFile = args.outSam,
    outFile = args.nameSortBam,
    name = True,
    threads = str(args.threads),
    memory = '2G',
    delete = True,
    path = 'samtools'
)
# Merge commands and run
alignSort = '%s && %s' %(alignCommand, sortCommand)
subprocess.check_output(alignSort, shell = True, stderr=subprocess.STDOUT)

###############################################################################
## Extract aligned pairs
###############################################################################
# extract pairs from alignments
alignMetrics, pairMetrics = alignedPair.extractPairs(
    inBam = args.nameSortBam,
    pairOut = args.outPairs,
    minMapQ = args.minMapQ,
    rmDup = args.rmDuplicates,
    rmConcord = args.rmConcordant,
    maxSize = args.maxSize
)
# Print alignment metrics
print '\nAlignment Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'total: %s' %(alignMetrics['total']),
    'unmapped: %s' %(alignMetrics['unmapped']),
    'poorly mapped: %s' %(alignMetrics['poormap']),
    'secondary alignment: %s' %(alignMetrics['secondary']),
    'singletons: %s' %(alignMetrics['singletons']),
    'multiple alignments: %s' %(alignMetrics['multiple']),
    'paired: %s' %(alignMetrics['pairs'])
)
# Print pair metrics
print '\nPair Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s' %( 
    'total: %s' %(pairMetrics['total']),
    'unique: %s' %(pairMetrics['unique']),
    'discordant: %s' %(pairMetrics['discord']),
    'unique discordant: %s' %(pairMetrics['discorduni']),
    'duplication rate: %.4f' %(1 - (pairMetrics['unique'] /
        float(pairMetrics['total']))),
    'concordant rate: %.4f' %(1 - (pairMetrics['discorduni'] /
        float(pairMetrics['unique'])))
)

##############################################################################
## Extract fragend pairs
##############################################################################
fragendMetrics = fragendPair.fragendPairs(
    pairIn = args.outPairs,
    fragendOut = args.outFrags,
    fasta = args.bwaFasta,
    maxDistance = args.maxDistance,
    resite = args.cutSite
)
# Print fragend metrics
print '\nFragend Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s' %(
    'total: %s' %(fragendMetrics['total']),
    'no fragend: %s' %(fragendMetrics['none']),
    'too distant: %s' %(fragendMetrics['distant']),
    'interchromosomal: %s' %(fragendMetrics['interchromosomal']),
    'intrachromosomal: %s' %(fragendMetrics['intrachromosomal'])
)
# Extract fragned and ligation metrics
fragendDist = fragendMetrics['fragDist']
ligationDist = fragendMetrics['ligDist']
fragendCount = len(fragendDist)
ligationCount = len(ligationDist)
# Print fragend metrics
if fragendCount > 0:
    fragendHist = numpy.histogram(
        fragendDist,
        bins = [0,100,200,300,400,500,600,700,800,900,1000,float('inf')]
    )
    fragendCumFreq = numpy.cumsum(fragendHist[0]) / float(fragendCount)
    # Print fragend distance data
    print '\nFragend Distance:'
    for bin, freq in itertools.izip(fragendHist[1][1:], fragendCumFreq):
        print '\t<%.0f: %.4f' %(
            bin,
            freq
        )
# Print ligation metrics
if ligationCount > 0:
    # Process ligation data
    print '\nIntrachromosomal Ligation Data:\n\t%s\n\t%s\n\t%s\n\t%s' %(
        'mean distance: %.0f' %(numpy.mean(ligationDist)),
        '25th percentile: %.0f' %(numpy.percentile(ligationDist,25)),
        '50th percentile: %.0f' %(numpy.percentile(ligationDist,50)),
        '75th percentile: %.0f' %(numpy.percentile(ligationDist,75))
    )
