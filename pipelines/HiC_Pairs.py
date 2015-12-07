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
from ngs_analysis.bam import samtools
from ngs_analysis.structure import alignedPair, fragendPair
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
print 'Parameters:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'minimum read length: ' + str(args.minLength),
    'mimimum mapping quality: ' + str(args.minMapQ),
    'maximum fragend distance: ' + str(args.maxDistance),
    'maximum distance of concordant pairs: ' + str(args.maxConcordant),
    'remove concordant pairs: ' + str(args.removeConcordant)
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
print 'Trim Metrics:\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
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
pairs, alignMetrics = alignedPair.extractPairs(
    inBam = args.nameSortBam,
    minMapQ = args.minMapQ
)
# Print alignment metrics
print 'Alignment Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'total: ' + str(alignMetrics['total']),
    'unmapped: ' + str(alignMetrics['unmapped']),
    'poorly mapped: ' + str(alignMetrics['poormap']),
    'secondary alignment: ' + str(alignMetrics['secondary']),
    'singletons: ' + str(alignMetrics['singletons']),
    'multiple alignments: ' + str(alignMetrics['multiple']),
    'paired: ' + str(alignMetrics['pairs'])
)
# Process duplicate and concordant reads
pairMetrics = alignedPair.processPairs(
    pairs = pairs,
    pairOut = args.outPairs,
    rmDup = True,
    rmConcord = args.removeConcordant,
    maxSize = args.maxDistance
)
# Print pair statistics
print 'Pair Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %( 
    'total: ' + str(pairMetrics['total']),
    'unique: ' + str(pairMetrics['unique']),
    'discordant: ' + str(pairMetrics['discord']),
    'unique discordant: ' + str(pairMetrics['discorduni']),
    'duplication rate: ' + pairMetrics['dupratio'],
    'concordant rate: ' + pairMetrics['conratio']
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
# Remove fragend and ligation distance data
fragendDist = fragendMetrics.pop('fragDist')
ligationDist = fragendMetrics.pop('ligDist')
# Print fragend metrics
print 'Fragend Data:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n' %(
    'total: ' + str(fragendMetrics['total']),
    'no fragend: ' + str(fragendMetrics['none']),
    'too distant: ' + str(fragendMetrics['distant']),
    'interchromosomal: ' + str(fragendMetrics['interchromosomal']),
    'intrachromosomal: ' + str(fragendMetrics['intrachromosomal'])
)
## Calclate and print fragend metrics if accepted fragends are present
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
