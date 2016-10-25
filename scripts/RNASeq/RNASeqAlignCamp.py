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
from general_python import slurm, docopt
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
if not os.path.isfile(args['<gtf>']):
    raise IOError('Could not find GTF file')
if not os.path.isfile(args['<rrna>']):
    raise IOError('Could not find rRNA list file')
# Convert numerical arguments
args['--threads'] = int(args['--threads'])
args['--forprob'] = float(args['--forprob'])
args['--minlength'] = int(args['--minlength'])
args['--trimqual'] = int(args['--trimqual'])
# Generate and store standard output directories
args['fastqDir'] = os.path.join(args['<outdir>'], 'trimFastq')
args['fastqSampleDir'] = os.path.join(args['fastqDir'], args['name'])
args['rsemTranDir'] = os.path.join(args['<outdir>'],  'rsemAlign')
args['rsemTranSampleDir'] = os.path.join(args['rsemTranDir'], args['name'])
args['tophatDir'] = os.path.join(args['<outdir>'], 'tophat2Align')
args['tophatSampleDir'] = os.path.join(args['tophatDir'], args['name'])
dirList = [args['fastqDir'], args['fastqSampleDir'], args['rsemTranDir'],
    args['rsemTranSampleDir'], args['tophatDir'], args['tophatSampleDir']]
# Create output directories
for directory in dirList:
    if not os.path.exists(directory):
        os.mkdir(directory)
# Extract path and module data and create object to submit commands
pmDict = slurm.parsePathModule(args['<paths>'])
slurmJobs = slurm.submitJobs()

###############################################################################
## Find and process fastq files
###############################################################################
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
# Check input read length
if len(args['read1']) < 1:
    raise IOError('Insufficient number of FASTQ files identified')
if len(args['read1']) == 1:
    args['read1'] = args['read1'][0]
    args['read2'] = args['read2'][0]
else:
    raise IOError('Need to develop multi read input')

###############################################################################
## Generate commands to trim fastq files
###############################################################################
# Generate file names
args['trimRead1'] = os.path.join(args['fastqSampleDir'], 
    args['name'] + '_trim_R1.fastq.gz')
args['trimRead2'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_trim_R2.fastq.gz')
args['trimLog'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_cutadapt.log')
# Generate commands to process single-end FASTQ files
trimReadsCommand = fastqTrim.cutadapt(
    read1In = args['read1'],
    read1Out = args['trimRead1'],
    read2In = args['read2'],
    read2Out = args['trimRead2'],
    path = pmDict[('cutadapt', 'path')],
    length = args['--minlength'],
    quality = args['--trimqual']
)
# Add trim command to path
trimReadsJobID = slurmJobs.add(
    command = trimReadsCommand,
    processors = 1,
    stdout = args['trimLog'],
    stderr = args['trimLog'],
    modules = pmDict[('cutadapt', 'modules')]
)

###############################################################################
## Generate commands to perform fastqc
###############################################################################
# Generate log file name
args['qcLog'] = os.path.join(args['fastqSampleDir'],
    args['name'] + '_fastqc.log')
# Create fastqc command
fastqcCommand = []
for fastq in (args['trimRead1'], args['trimRead2']):
    if fastq is None:
        continue
    fastqcCommand.append(
        fastqQC.fastQC(
            inFile = fastq, 
            outDir = args['fastqSampleDir'],
            path = pmDict[('fastqc', 'path')]
        )
    )
fastqcCommand = '; '.join(fastqcCommand)
# Submit FASTQC commands
fastqcJobID = slurmJobs.add(
    command = fastqcCommand,
    stdout = args['qcLog'],
    stderr = args['qcLog'],
    depend = [trimReadsJobID],
    modules = pmDict[('fastqc', 'modules')]
)

###############################################################################
## Generate RSEM Transcript alignment command
###############################################################################
# Generate file names
args['rsemTranPrefix'] = os.path.join(args['rsemTranSampleDir'], args['name'])
args['rsemTranLog'] = args['rsemTranPrefix'] + '_rsem.log'
args['rsemTranBam'] = args['rsemTranPrefix'] + '.transcript.bam'
args['rsemTranGenomeBam'] = args['rsemTranPrefix'] + '.genome.bam'
# Generate rsem transcript alignment command
rsemAlignCommand = fastqAlign.rsemBowtie2Align(
    read1 = args['trimRead1'],
    read2 = args['trimRead2'],
    index = args['<rsemtranindex>'],
    outPrefix = args['rsemTranPrefix'],
    rsemPath = pmDict[('rsem', 'path')],
    threads = args['--threads'],
    forProb = args['--forprob'],
    genomeBam = False
)
# Submit rsem transcript command
rsemAlignJobID = slurmJobs.add(
    command = rsemAlignCommand,
    processors = args['--threads'],
    stdout = args['rsemTranLog'],
    stderr = args['rsemTranLog'],
    depend = [trimReadsJobID],
    modules = pmDict[('rsem', 'modules')]
)

###############################################################################
## Generate Tophat2 Alignment Command
###############################################################################
# Generate file names
args['tophatLog'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_tophat2.log')
# Generate tophat command
tophat2AlignCommand = fastqAlign.tophat2Align(
    read1 = args['trimRead1'],
    read2 = args['trimRead2'],
    genomeIndex = args['<bwt2genindex>'],
    transcriptIndex = args['<bwt2tranindex>'],
    outDir = args['tophatSampleDir'],
    path = pmDict[('tophat2', 'path')],
    forProb = args['--forprob'],
    threads = args['--threads'],
    readGroup = '1',
    sampleName = args['name']
)
# Submit tophat2 commands
tophat2AlignJobID = slurmJobs.add(
    command = tophat2AlignCommand,
    processors = args['--threads'],
    stdout = args['tophatLog'],
    stderr = args['tophatLog'],
    depend = [trimReadsJobID],
    modules = pmDict[('tophat2', 'modules')]
)

###############################################################################
## Mark duplicates in Tophat output
###############################################################################
args['tophatBam'] = os.path.join(args['tophatSampleDir'], 'accepted_hits.bam')
args['mdupBam'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_tophat2_mdup.bam')
args['mdupLog1'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_mdup.log1')
args['mdupLog2'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_mdup.log2')
# Create command to mark duplicates
markDupCommand = picard.markDuplicates(
    inBam = args['tophatBam'],
    outBam = args['mdupBam'],
    logFile = args['mdupLog1'],
    removeDuplicates = False,
    picardPath = pmDict[('picard', 'path')],
    delete = True,
    memory = 10
)
# Submit mark duplicates job ID
markDupJobID = slurmJobs.add(
    command = markDupCommand,
    stdout = args['mdupLog2'],
    stderr = args['mdupLog2'],
    depend = [tophat2AlignJobID],
    memory = 12,
    modules = pmDict[('picard', 'modules')]
)

###############################################################################
## Submit RNASeqC command
###############################################################################
args['seqcLog'] = os.path.join(args['tophatSampleDir'],
    args['name'] + '_rnaseqc.log')
# Create command to peform RNASeqC
seqcCommand = bamQC.RNASeqC(
    inBam = args['mdupBam'],
    fasta = args['<bwt2genindex>'] + '.fa',
    gtf = args['<gtf>'],
    rRNA = args['<rrna>'],
    outDir = args['tophatSampleDir'],
    outPrefix = args['name'],
    seqcPath = pmDict[('rnaseqc', 'path')],
    singleEnd = args['--singleend'],
    memory = 10
)
# Submit RNASeqC command
seqcJobID = slurmJobs.add(
    command = seqcCommand,
    stdout = args['seqcLog'],
    stderr = args['seqcLog'],
    depend = [markDupJobID],
    memory = 12,
    modules = pmDict[('rnaseqc', 'modules')]
)

#################################################################################
## Submit commands
#################################################################################
# Submit jobs and print command
slurmJobs.submit(verbose = True)
