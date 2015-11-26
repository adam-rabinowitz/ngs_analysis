#!/farm/babs/redhat6/software/python-2.7.3/bin/python
''' A script to quality trim and align RNA-Seq data. Script aligns the read data
to both the transcriptome and spike-in sequences using RSEM. The script also
peforms a Tophat2 alignment to the genome. To determine the quality of the
dataset an RNASeqC analysis is performed on the results of the Tophat2 output.
Script takes X arguments.
1)  sampleData - A string listing the sample name, library name and read
    group of the data. Elements should be delimited by the ':' character. The
    sample name will be used as a prefix of the results files. The library name
    will be used as the prefix with which to identify the FASTQ files for the
    analysis.
2)  inDir - A comma seperated list of the directories containing the FASTQ
    files.
3)  outDir - Path to the output directory. If they do not already exist, folders
    titled 'fastqData', 'rsemTranscriptData', 'rsemSpikeData' and 'tophatData
    will be created within the output directory.
4)  rsemIndex - Prefix of Bowtie2 index of transcriptome. Created from
    genomic fasta file and gtf file using the RSEM command
    'rsem-prepare_reference'.
5)  bowtie2GenomeIndex - Prefix of Bowtie2 index of genome fasta file. Created
    using Bowtie2 'bowtie2-build' command. A genome fasta file of the same
    prefix should pe present. A sequence dictionary and index of the fasta
    file should also be present.
6)  bowtie2TranscriptIndex - Prefix''' 
################################################################################
## Load modules and extract command line arguments
################################################################################
# Import functions from external modules
import subprocess
import sys
import argparse
import re
import os
# Import custom modules
from ngs_analysis.fastq import fastqTrim, fastqQC, fastqAlign, fastqExtract, fastqFind
from ngs_analysis.bam import samtools, picard, bamQC
from ngs_analysis.system import moab
import alignment_functions
# Print command
print '%s\n' %(' '.join(sys.argv))
# Process command line options
parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
)
# Define positional arguments
parser.add_argument(
    'sampleData',
    help = 'Comma seperated list of Name and FASTQ prefix of sample',
    type = str
)
parser.add_argument(
    'inDir',
    help = 'Comma seperated list of FASTQ directories',
    type = str
)
parser.add_argument(
    'outDir',
    help = 'Path to output directory',
    type = str
)
parser.add_argument(
    'rsemTranIndex',
    help = 'Suffix of RSEM transcript index',
    type = str
)
parser.add_argument(
    'bowtie2GenomeIndex',
    help = 'Suffix of Bowtie2 genomic index',
    type = str
)
parser.add_argument(
    'bowtie2TranscriptIndex',
    help = 'Suffix of Tophat2 transcript index',
    type = str
)
parser.add_argument(
    'gtfFile',
    help = 'GTF file',
    type = str
)
parser.add_argument(
    'rRNAList',
    help = 'List file containing intervals covered by rRNA transcripts',
    type = str
)
# Define optional arguments
parser.add_argument('-i', '--rsemSpikeIndex', help = 'Suffix of RSEM'+\
    ' spike-in index', default = '', type = str)
parser.add_argument('-f', '--forwardProbability', help = 'Probability of'+\
    ' first read deriving from forward strand', default = '0.5', type = str)
parser.add_argument('-q', '--trimQuality', help = 'Minimum Base quality for'+\
    ' trimming', default = '20', type = str)
parser.add_argument(
    '-l', '--minLength', help = 'Minimum read length post trimming',
    default = '25', type = str)
parser.add_argument('-m', '--metrics_script', help = 'Script to calculate '+\
    'tophat2 mate distance metrics', default = '/farm/home/rabino01/'+\
    'python/calculate_insert_metrics.py', type = str)
parser.add_argument('-t', '--tophat2', help = 'Path to Tophat2 executable',
    default = 'tophat', type = str)
parser.add_argument('-b', '--bowtie2', help = 'Path to Bowtie2 executable',
    default = 'bowtie2', type = str)
parser.add_argument('-r', '--rsem', help = 'Path to rsem-calculate-expression',
    default = 'rsem-calculate-expression',  type = str)
parser.add_argument('-s', '--samtools', help = 'Path to samtools executable',
    default = 'samtools',type = str)
parser.add_argument('-n', '--rnaseqc', help = 'RNASeqC jar file',
    default = 'RNASeQC.jar', type = str)
parser.add_argument('-p', '--picard', help = 'Picard jar file',
    default = 'picard.jar', type = str)
parser.add_argument('-a', '--cutadapt', help = 'Path to cutadapt executable',
    default = 'cutadapt', type = str)
parser.add_argument('-c', '--fastqc', help = 'FastQC exectuable',
    default = 'fastqc', type = str)
parser.add_argument('-j', '--java', help = 'Java Path',
    default = 'java', type = str)
# Extract arguments
args = parser.parse_args()

################################################################################
## Process command line arguments and create output directories
################################################################################
# Extract sample data
args.name, args.prefix = args.sampleData.split(',')
# Extract fastq files
args.read1, args.read2 = fastqFind.findIlluminaFastq(
    prefix = args.prefix,
    dirList = args.inDir.split(','),
    pair = True
)
# Check number of FASTQ files
if len(args.read1) != len(args.read2) or len(args.read1) < 1:
    raise IOError('Unequal or insufficient number of FASTQ files identified')
# Process strand arguments
if float(args.forwardProbability) < 0 or float(args.forwardProbability) > 1:
    raise ValueError("Forward probability should be 1 >= and >= 0")
# Check python script to calculate alignment metrics exists
if not os.path.isfile(args.metrics_script):
    raise IOError('Python script to calculate insert metrics does not exist')
# Generate and store standard output directories
args.fastqDir = args.outDir + 'fastq/'
args.fastqSampleDir = args.fastqDir + args.name + '/'
args.rsemTranDir = args.outDir + 'rsemTran/'
args.rsemTranSampleDir = args.rsemTranDir + args.name + '/'
args.tophatDir = args.outDir + 'tophat/'
args.tophatSampleDir = args.tophatDir + args.name + '/'
dirList = [args.fastqDir, args.fastqSampleDir, args.rsemTranDir,
    args.rsemTranSampleDir, args.tophatDir, args.tophatSampleDir]
# Generate and store spike-in output directories if argument supplied
if args.rsemSpikeIndex:
    args.rsemSpikeDir = args.outDir + 'rsemSpike/'
    args.rsemSpikeSampleDir = args.rsemSpikeDir + args.name + '/'
    dirList.extend([args.rsemSpikeDir, args.rsemSpikeSampleDir])
# Create output directories
for directory in dirList:
    if not os.path.exists(directory):
        os.mkdir(directory)

#################################################################################
### Generate commands to concatenate and trim FASTQ files
#################################################################################
# Generate fastq file names
args.concatRead1 = args.fastqSampleDir + args.name + '_concat_R1.fastq.gz'
args.concatRead2 = args.fastqSampleDir + args.name + '_concat_R2.fastq.gz'
args.trimRead1 = args.fastqSampleDir + args.name + '_trim_R1.fastq.gz'
args.trimRead2 = args.fastqSampleDir + args.name + '_trim_R2.fastq.gz'
args.trimLog = args.fastqSampleDir + args.name + '_trim.log'
# Generate commands to process mutliple FASTQ files
if len(args.read1) > 1:
    # Create commands to concatenate files
    concatRead1Command = 'zcat %s | gzip > %s' %(
        ' '.join(args.read1),
        args.concatRead1
    )
    concatRead2Command = 'zcat %s | gzip > %s' %(
        ' '.join(args.read2),
        args.concatRead2
    )
    # Create command to trim FASTQ files
    trim = fastqTrim.cutadaptTrimPaired(
        read1In = args.concatRead1,
        read2In = args.concatRead2,
        read1Out = args.trimRead1,
        read2Out = args.trimRead2,
        path = args.cutadapt
    )
    trimReadsCommand = '%s; rm %s %s' %(
        trim,
        args.concatRead1,
        args.concatRead2
    )
# Generate commands to process single FASTQ files
else:
    concatRead1Command = ''
    concatRead2Command = ''
    trimReadsCommand = fastqTrim.cutadaptTrimPaired(
        read1In = args.read1[0],
        read2In = args.read2[0],
        read1Out = args.trimRead1,
        read2Out = args.trimRead2,
        path = args.cutadapt
    )

################################################################################
## Perform FastQC on trimmed reads
################################################################################
# Generate commands
read1Fastqc = fastqQC.fastQC(
    inFile = args.trimRead1, 
    outDir = args.fastqSampleDir,
    path = args.fastqc
)
read2Fastqc = fastqQC.fastQC(
    inFile = args.trimRead2,
    outDir = args.fastqSampleDir,
    path = args.fastqc
)
# Concatenate commands
fastqcCommand = '%s && %s' %(
    read1Fastqc,
    read2Fastqc
)

################################################################################
## Generate RSEM Transcript alignment command
################################################################################
# Generate file names
args.rsemTranPrefix = args.rsemTranSampleDir + args.name
args.rsemTranLog = args.rsemTranPrefix + '_RSEM.log'
args.rsemTranBam = args.rsemTranPrefix + '.transcript.bam'
args.rsemTranGenomeBam = args.rsemTranPrefix + '.genome.bam'
# Generate rsem transcript alignment command
rsemTranAlign = fastqAlign.rsemBowtie2Align(
    read1 = args.trimRead1,
    read2 = args.trimRead2,
    index = args.rsemTranIndex,
    outPrefix = args.rsemTranPrefix,
    rsemPath = args.rsem,
    bowtie2Path = re.sub('[^/]*$', '', args.bowtie2),
    threads = 4,
    forProb = args.forwardProbability,
    genomeBam = True
)
rsemTranCommand = '%s && rm %s %s' %(
    rsemTranAlign,
    args.rsemTranBam,
    args.rsemTranGenomeBam
)

################################################################################
## Generate RSEM Spike-In Alignment command
################################################################################
# Generate file names
if args.rsemSpikeIndex:
    args.rsemSpikePrefix = args.rsemSpikeSampleDir + args.name
    args.rsemSpikeLog = args.rsemSpikePrefix + '_RSEM.log'
    args.rsemSpikeBam = args.rsemSpikePrefix + '.transcript.bam'
    # Generate rsem transcript alignment command
    rsemSpikeAlign = fastqAlign.rsemBowtie2Align(
        read1 = args.trimRead1,
        read2 = args.trimRead2,
        index = args.rsemSpikeIndex,
        outPrefix = args.rsemSpikePrefix,
        rsemPath = args.rsem,
        bowtie2Path = re.sub('[^/]*$', '', args.bowtie2),
        threads = 4,
        forProb = args.forwardProbability,
        genomeBam = False
    )
    rsemSpikeCommand = '%s && rm %s' %(
        rsemSpikeAlign,
        args.rsemSpikeBam
    )

################################################################################
## Generate Tophat2 Alignment Command
################################################################################
# Generate file names
args.sampleRead1 = args.tophatSampleDir + args.name + '_sample_R1.fastq.gz'
args.sampleRead2 = args.tophatSampleDir + args.name + '_sample_R2.fastq.gz'
args.tempSam = args.tophatSampleDir + args.name + '_temp.sam'
args.tophatLog = args.tophatSampleDir + args.name + '_Tophat2.log'
# Extract sample command
sampleCommand = fastqExtract.extractRandom(
    read1In = args.trimRead1,
    read2In = args.trimRead2,
    read1Out = args.sampleRead1,
    read2Out = args.sampleRead2,
    number = 250000
)
# Extract alignment command
alignSampleCommand = fastqAlign.bowtie2Align(
    read1 = args.sampleRead1,
    read2 = args.sampleRead2,
    index = args.bowtie2TranscriptIndex,
    path = args.bowtie2,
    outSam = args.tempSam,
    threads = 4
)
# Calculate insert metrics
metricsCommand = 'read m s <<< $(python %s %s)' %(
    args.metrics_script,
    args.tempSam
)
# Generate tophat2 command
tophat2Command = fastqAlign.tophat2Align(
    read1 = args.trimRead1,
    read2 = args.trimRead2,
    genomeIndex = args.bowtie2GenomeIndex,
    transcriptIndex = args.bowtie2TranscriptIndex,
    outDir = args.tophatSampleDir,
    path = args.tophat2,
    forProb = args.forwardProbability,
    threads = 4,
    mateDist = '$m',
    mateSD = '$s',
    readGroup = args.name,
    sampleName = args.name
)
# Combine commands
tophat2Combined = '%s && %s && %s && %s && rm %s %s' %(
    sampleCommand,
    alignSampleCommand,
    metricsCommand,
    tophat2Command,
    args.sampleRead1,
    args.sampleRead2
)

################################################################################
## Perform RNASeQC analysis
################################################################################
args.tophatBam = args.tophatSampleDir + 'accepted_hits.bam'
args.mdupBam = args.tophatSampleDir + args.name + '_tophat_mdup.bam'
args.mdupLog = args.tophatSampleDir + args.name + '_mdup.log'
args.seqcLog = args.tophatSampleDir + args.name + '_rnaseqc.log'
# Create command to mark duplicates
markDupCommand = picard.markDuplicates(
    inBam = args.tophatBam,
    outBam = args.mdupBam,
    logFile = args.mdupLog,
    removeDuplicates = False,
    picardPath = args.picard,
    javaPath = args.java,
    delete = True
)
# Create command to peform RNASeqC
seqcCommand = bamQC.RNASeqC(
    inBam = args.mdupBam,
    fasta = args.bowtie2GenomeIndex + '.fa',
    gtf = args.gtfFile,
    rRNA = args.rRNAList,
    outDir = args.tophatSampleDir,
    outPrefix = args.name,
    seqcPath = args.rnaseqc,
    javaPath = args.java
)
# Combine mark, index and RNASeqC commands
seqcComboCommand = '%s && %s' %(
    markDupCommand,
    seqcCommand
)

################################################################################
# Submit commands
################################################################################
#Submit concatenation commands if required
if concatRead1Command:
    concatRead1JobID = moab.submitjob(
        concatRead1Command
    )
    concatRead2JobID = moab.submitjob(
        concatRead2Command
    )
    print 'Concatenate read1 command: %s' %(concatRead1Command)
    print 'Concatenate read1 job ID: %s\n' %(concatRead1JobID)
    print 'Concatenate read2 command: %s' %(concatRead2Command)
    print 'Concatenate read2 job ID: %s\n' %(concatRead2JobID)
    concatReadsJobID = [concatRead1JobID, concatRead2JobID]
else:
    concatReadsJobID = None
# Submit trimming commands
trimReadsJobID = moab.submitjob(
    command = trimReadsCommand,
    stdout = args.trimLog,
    stderr = args.trimLog,
    dependency = concatReadsJobID
)
print 'Trim command: %s' %(trimReadsCommand)
print 'Trim Job ID: %s\n' %(trimReadsJobID)
# Submit FASTQC commands
fastqcJobID = moab.submitjob(
    command = fastqcCommand,
    processor = 1,
    dependency = trimReadsJobID
)
print 'FastQC Command: %s' %(fastqcCommand)
print 'FastQC Job ID: %s\n' %(fastqcJobID)
# Submit rsem transcript command
rsemTranJobID = moab.submitjob(
    command = rsemTranCommand,
    processor = 4,
    stdout = args.rsemTranLog,
    stderr = args.rsemTranLog,
    dependency = trimReadsJobID
)
print 'RSEM Transcript Command: %s' %(rsemTranCommand)
print 'RSEM Transcript Job ID: %s\n' %(rsemTranJobID)
# Submit rsem spike command
if args.rsemSpikeIndex:
    rsemSpikeJobID = moab.submitjob(
        command = rsemSpikeCommand,
        processor = 4,
        stdout = args.rsemSpikeLog,
        stderr = args.rsemSpikeLog,
        dependency = trimReadsJobID
    )
    print 'RSEM Spike Command: %s' %(rsemSpikeCommand)
    print 'RSEM Spike Job ID: %s\n' %(rsemSpikeJobID)
# Submit tophat2 commands
tophatJobID = moab.submitjob(
    command = tophat2Combined,
    processor = 4,
    stdout = args.tophatLog,
    stderr = args.tophatLog,
    dependency = trimReadsJobID
)
print 'Tophat Command: %s' %(tophat2Combined)
print 'Tophat Job ID: %s\n' %(tophatJobID)
# Submit RNASeqC command
seqcJobID = moab.submitjob(
    command = seqcComboCommand,
    processor = 1,
    stdout = args.seqcLog,
    stderr = args.seqcLog,
    dependency = tophatJobID
)
print 'RNASeqC Command: %s' %(seqcComboCommand)
print 'RNASeqC JobID: %s' %(seqcJobID)
