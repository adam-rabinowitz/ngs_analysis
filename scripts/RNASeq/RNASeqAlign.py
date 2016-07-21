'''RNASeqAlign.py

Usage:
    
    RNASeqAlign.py <sampledata> <indir> <outdir> <rsemtranindex>
        <bwt2genindex> <bwt2tranindex> <gtf> <rrna> <paths>
        [--rsemspikeindex=<rsi>] [--forprob=<fp>]
        [--trimqual=<tq>] [--minlength=<ml>] [--threads=<th>]
        [--singleend]
    
Options:
    
    --rsemspikeindex=<rsi>  RSEM index for spike-in transcripts
    --forprob=<fp>          Probability read1 from sense strand [default: 0.5]
    --trimqual=<tq>         Base quality for trimming [default: 20]
    --minlength=<ml>        Minimum length of reads [default: 25]
    --threads=<th>          Number of threads to use [default: 4]
    --singleend             Samples are single end

'''
###############################################################################
## Load modules and extract command line arguments
###############################################################################
# Import functions from external modules
import subprocess
import sys
import argparse
import re
import os
# Import custom modules
from ngs_python.fastq import fastqTrim, fastqQC, fastqAlign, fastqFind
from ngs_python.bam import samtools, picard, bamQC
from general_python import moab, docopt, toolbox
# Print command
print '%s\n' %(' '.join(sys.argv))

###############################################################################
## Process command line arguments and create output directories
###############################################################################
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Extract sample prefix and name = args['<sampledata>'].split(',')
args['prefix'], args['name'] = args['<sampledata>'].split(',')
# Check supplied files
toolbox.check_var(args['<gtf>'], 'file')
toolbox.check_var(args['<rrna>'], 'file')
# Extract fastq files and check
if args['--singleend']:
    args['read1']  = fastqFind.findFastq(
        prefix = args['prefix'],
        dirList = args['<indir>'].split(','),
        pair = False
    )
else:
    args['read1'], args['read2'] = fastqFind.findFastq(
        prefix = args['prefix'],
        dirList = args['<indir>'].split(','),
        pair = True
    )
    if len(args['read1']) != len(args['read2']):
        raise IOError('Unequal number of FASTQ files identified')
if len(args['read1']) < 1:
    raise IOError('Insufficient number of FASTQ files identified')
# Convert numerical arguments
args['--threads'] = int(args['--threads'])
args['--forprob'] = float(args['--forprob'])
args['--minlength'] = int(args['--minlength'])
args['--trimqual'] = int(args['--trimqual'])
# Generate and store standard output directories
args['fastqDir'] = os.path.join(args['<outdir>'], 'fastq')
args['fastqSampleDir'] = os.path.join(args['fastqDir'], args['name'])
args['rsemTranDir'] = os.path.join(args['<outdir>'],  'rsemTran')
args['rsemTranSampleDir'] = os.path.join(args['rsemTranDir'], args['name'])
args['tophatDir'] = os.path.join(args['<outdir>'], 'tophat')
args['tophatSampleDir'] = os.path.join(args['tophatDir'], args['name'])
dirList = [args['fastqDir'], args['fastqSampleDir'], args['rsemTranDir'],
    args['rsemTranSampleDir'], args['tophatDir'], args['tophatSampleDir']]
# Generate and store spike-in output directories if argument supplied
if args['--rsemspikeindex']:
    args['rsemSpikeDir'] = os.path.join(args['<outdir>'] + 'rsemSpike')
    args['rsemSpikeSampleDir'] = os.path.join(args['rsemSpikeDir'],
        args['name'])
    dirList.extend([args['rsemSpikeDir'], args['rsemSpikeSampleDir']])
# Create output directories
for directory in dirList:
    if not os.path.exists(directory):
        os.mkdir(directory)
# Extract path data and check
paths = toolbox.fileDict(args['<paths>'], sep = '\t')
for program in paths:
    toolbox.check_var(paths[program], 'exc')
# Generate object to store and submit commands
moabJobs = moab.moabJobs()

#################################################################################
### Generate commands to trim and QC fastq files
#################################################################################
# Generate fastq file names
args['trimRead1'] = os.path.join(args['fastqSampleDir'], 
    args['name'] + '_trim_R1.fastq.gz')
args['trimRead2'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_trim_R2.fastq.gz')
args['trimLog'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_trim.log')
args['qcLog'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_qc.log')
# Generate commands to process single-end FASTQ files
if args['--singleend']:
    # Create trim command
    trimReadsCommand = fastqTrim.cutadaptTrim(
        readIn = args['read1'],
        readOut = args['trimRead1'],
        path = paths['cutadapt'],
        length = args['--minlength'],
        quality = args['--trimqual']
    )
    # Create QC command
    fastqcCommand = fastqQC.fastQC(
        inFile = args['trimRead1'], 
        outDir = args['fastqSampleDir'],
        path = paths['fastqc']
    )
    # Empty read2 file
    args['trimRead2'] = None
# Generate commands to process paired-end FASTQ files
else:
    # Create trim command
    trimReadsCommand = fastqTrim.cutadaptTrimPaired(
        read1In = args['read1'],
        read2In = args['read2'],
        read1Out = args['trimRead1'],
        read2Out = args['trimRead2'],
        path = paths['cutadapt'],
        length = args['--minlength'],
        quality = args['--trimqual']
    )
    # Generate QC command
    read1Fastqc = fastqQC.fastQC(
        inFile = args['trimRead1'], 
        outDir = args['fastqSampleDir'],
        path = paths['fastqc']
    )
    read2Fastqc = fastqQC.fastQC(
        inFile = args['trimRead2'],
        outDir = args['fastqSampleDir'],
        path = paths['fastqc']
    )
    fastqcCommand = '%s && %s' %(
        read1Fastqc,
        read2Fastqc
    )
# Submit trimming commands
trimReadsJobID = moabJobs.add(
    command = trimReadsCommand,
    stdout = args['trimLog'],
    stderr = args['trimLog']
)
# Submit FASTQC commands
fastqcJobID = moabJobs.add(
    command = fastqcCommand,
    stdout = args['qcLog'],
    stderr = args['qcLog'],
    dependency = [trimReadsJobID]
)

################################################################################
## Generate RSEM Transcript alignment command
################################################################################
# Generate file names
args['rsemTranPrefix'] = os.path.join(args['rsemTranSampleDir'], args['name'])
args['rsemTranLog'] = args['rsemTranPrefix'] + '_RSEM.log'
args['rsemTranBam'] = args['rsemTranPrefix'] + '.transcript.bam'
args['rsemTranGenomeBam'] = args['rsemTranPrefix'] + '.genome.bam'
# Generate rsem transcript alignment command
rsemAlignCommand = fastqAlign.rsemBowtie2Align(
    read1 = args['trimRead1'],
    read2 = args['trimRead2'],
    index = args['<rsemtranindex>'],
    outPrefix = args['rsemTranPrefix'],
    rsemPath = paths['rsem'],
    bowtie2Path = re.sub('[^/]*$', '', paths['bowtie2']),
    threads = args['--threads'],
    forProb = args['--forprob'],
    genomeBam = True
)
# Submit rsem transcript command
rsemAlignJobID = moabJobs.add(
    command = rsemAlignCommand,
    processors = args['--threads'],
    stdout = args['rsemTranLog'],
    stderr = args['rsemTranLog'],
    dependency = [trimReadsJobID]
)

################################################################################
## Generate RSEM Spike-In Alignment command
################################################################################
# Generate file names
if args['--rsemspikeindex']:
    args['rsemSpikePrefix'] = os.path.join(args['rsemSpikeSampleDir'],
        args['name'])
    args['rsemSpikeLog'] = args['rsemSpikePrefix'] + '_RSEM.log'
    args['rsemSpikeBam'] = args['rsemSpikePrefix'] + '.transcript.bam'
    # Generate rsem transcript alignment command
    rsemSpikeAlign = fastqAlign.rsemBowtie2Align(
        read1 = args['trimRead1'],
        read2 = args['trimRead2'],
        index = args['--rsemspikeindex'],
        outPrefix = args['rsemSpikePrefix'],
        rsemPath = paths['rsem'],
        bowtie2Path = re.sub('[^/]*$', '', paths['bowtie2']),
        threads = args['--threads'],
        forProb = args['--forprob'],
        genomeBam = False
    )
    rsemSpikeCommand = '%s && rm %s' %(
        rsemSpikeAlign,
        args['rsemSpikeBam']
    )
    # Submit command
    rsemSpikeJobID = moabJobs.add(
        command = rsemSpikeCommand,
        processors = args['--threads'],
        stdout = args['rsemSpikeLog'],
        stderr = args['rsemSpikeLog'],
        dependency = [trimReadsJobID]
    )

###############################################################################
## Generate Tophat2 Alignment Command
###############################################################################
# Generate file names
args['tempSam'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_sample.sam')
args['nsortBam'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_sample_nsort.bam')
args['tophatLog'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_Tophat2.log')
# Generate command for single-end samples
if args['--singleend']:
    # Generate command to perform tophat2 alignment
    tophat2Complete = fastqAlign.tophat2Align(
        read1 = args['trimRead1'],
        genomeIndex = args['<bwt2genindex>'],
        transcriptIndex = args['<bwt2tranindex>'],
        outDir = args['tophatSampleDir'],
        path = paths['tophat2'],
        forProb = args['--forprob'],
        threads = args['--threads'],
        readGroup = '1',
        sampleName = args['name']
    )
# Generate command for paired-end samples
else:
    # Generate command to align first 1e6 reads to transcriptome
    alignSampleCommand = fastqAlign.bowtie2Align(
        read1 = args['trimRead1'],
        read2 = args['trimRead2'],
        index = args['<bwt2tranindex>'],
        bowtie2Path = paths['bowtie2'],
        samtoolsPath = paths['samtools'],
        outFile = args['nsortBam'],
        threads = args['--threads'],
        upto = 1000000
    )
    # Calculate insert metrics
    metricsCommand = 'read m s <<< $(%s %s %s)' %(
        paths['python'],
        paths['distance'],
        args['nsortBam']
    )
    # Generate tophat2 command
    tophat2Command = fastqAlign.tophat2Align(
        read1 = args['trimRead1'],
        read2 = args['trimRead2'],
        genomeIndex = args['<bwt2genindex>'],
        transcriptIndex = args['<bwt2tranindex>'],
        outDir = args['tophatSampleDir'],
        path = paths['tophat2'],
        forProb = args['--forprob'],
        threads = args['--threads'],
        mateDist = '$m',
        mateSD = '$s',
        readGroup = '1',
        sampleName = args['name']
    )
    # Combine commands
    tophat2Complete = '%s && %s && %s' %(
        alignSampleCommand,
        metricsCommand,
        tophat2Command
    )
# Submit tophat2 commands
tophatJobID = moabJobs.add(
    command = tophat2Complete,
    processors = args['--threads'],
    stdout = args['tophatLog'],
    stderr = args['tophatLog'],
    dependency = [trimReadsJobID]
)

################################################################################
## Perform RNASeQC analysis
################################################################################
args['tophatBam'] = os.path.join(args['tophatSampleDir'], 'accepted_hits.bam')
args['mdupBam'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_tophat_mdup.bam')
args['mdupLog'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_mdup.log')
args['seqcLog'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_rnaseqc.log')
# Create command to mark duplicates
markDupCommand = picard.markDuplicates(
    inBam = args['tophatBam'],
    outBam = args['mdupBam'],
    logFile = args['mdupLog'],
    removeDuplicates = False,
    picardPath = paths['picard'],
    javaPath = paths['java'],
    delete = True
)
# Create command to peform RNASeqC
seqcCommand = bamQC.RNASeqC(
    inBam = args['mdupBam'],
    fasta = args['<bwt2genindex>'] + '.fa',
    gtf = args['<gtf>'],
    rRNA = args['<rrna>'],
    outDir = args['tophatSampleDir'],
    outPrefix = args['name'],
    seqcPath = paths['rnaseqc'],
    javaPath = paths['java'],
    singleEnd = args['--singleend']
)
# Combine mark, index and RNASeqC commands
seqcComboCommand = '%s && %s' %(
    markDupCommand,
    seqcCommand
)
# Submit RNASeqC command
seqcJobID = moabJobs.add(
    command = seqcComboCommand,
    processors = 1,
    stdout = args['seqcLog'],
    stderr = args['seqcLog'],
    dependency = [tophatJobID]
)

################################################################################
# Submit commands
################################################################################
# Submit jobs and print command
moabJobs.submit(verbose = True)
