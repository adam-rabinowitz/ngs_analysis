'''ATACSeqAlignBed.py

Usage:
    
    ATACSeqAlignBed.py <sampledata> <indir> <outdir> <bwt2index>
        <chrfile> <paths> [--minMapQ=<minMapQ>] [--interval=<interval>]
        [--threads=<threads>] [--rmDup]
    
Options:
    
    --minMapQ=<minMapQ>    Minimum read mapping quality [default: 20]
    --interval=<interval>  Read start interval size [default: 50]
    --threads=<threads>    Number of threads to use [default: 4]
    --rmDup                Flag to ignore duplicate alignments

'''


# Import required modules
import os
import sys
from ngs_python.fastq import fastqFind, fastqTrim, fastqQC, fastqAlign
from ngs_python.bam import picard
from ngs_python.bed import bedtools
from general_python import docopt, moab, toolbox
# Print sommand and extract argument
print '%s\n' %(' '.join(sys.argv))
args = docopt.docopt(__doc__,version = 'v1')
# Extract sample prefix and name
args['prefix'], args['name'] = args['<sampledata>'].split(',')
# Extract path data and check they are executable:
paths = toolbox.fileDict(args['<paths>'], sep = '\t')
for program in paths:
    toolbox.checkArg(paths[program], 'exc')
# Conver numric arguments
args['--minMapQ'] = int(args['--minMapQ'])
args['--threads'] = int(args['--threads'])
args['--interval'] = int(args['--interval'])
# Extract fastq files and check
args['read1'], args['read2'] = fastqFind.findFastq(
    prefix = args['prefix'],
    dirList = args['<indir>'].split(','),
    pair = True)
if len(args['read1']) != len(args['read2']):
    raise IOError('Unequal number of FASTQ files identified')
if len(args['read1']) < 1:
    raise IOError('Insufficient number of FASTQ files identified')
# Check output directories
if not os.path.isdir(args['<outdir>']):
    raise IOError('Output directory not found')
# Create directory names
args['fsqDir'] = os.path.join(args['<outdir>'], 'trimFastq')
args['bamDir'] = os.path.join(args['<outdir>'], 'bamFiles')
args['bedDir'] = os.path.join(args['<outdir>'], 'bedFiles')
args['smpFsq'] = os.path.join(args['fsqDir'], args['name'])
args['smpBam'] = os.path.join(args['bamDir'], args['name'])
args['smpBed'] = os.path.join(args['bedDir'], args['name'])
# Create directories
for directory in [args['fsqDir'], args['bamDir'], args['bedDir'],
    args['smpFsq'], args['smpBam'], args['smpBed']]:
    if not os.path.exists(directory):
        os.mkdir(directory)
# Create output files for FASTQ processing
args['fsqPrf'] = os.path.join(args['smpFsq'], args['name'])
args['trmRd1'] = args['fsqPrf'] + '.trim.R1.fastq.gz'
args['trmRd2'] = args['fsqPrf'] + '.trim.R2.fastq.gz'
args['trmLog'] = args['fsqPrf'] + '.trim.log'
args['fqcLog'] = args['fsqPrf'] + '.fastqc.log'
# Generate output files for BAM processing
args['bamPrf'] = os.path.join(args['smpBam'], args['name'])
args['alnBam'] = args['bamPrf'] + '.initial.bam'
args['alnLog'] = args['bamPrf'] + '.align.log'
args['mdpBam'] = args['bamPrf'] + '.bam'
args['mdpLog1'] = args['bamPrf'] + '.mdup.1.log'
args['mdpLog2'] = args['bamPrf'] + '.mdup.2.log'
args['insPrf'] = args['bamPrf'] + '.insert'
args['insLog'] = args['bamPrf'] + '.insert.log'
# Generate output files for BED processing
args['bedPrf'] = os.path.join(args['smpBed'], args['name'])
args['intBed'] = args['bedPrf'] + '.initial.bed'
args['intLog'] = args['bedPrf'] + '.bam2bed.log'
args['srtBed'] = args['bedPrf'] + '.bed'
args['bedGrp'] = args['bedPrf'] + '.bedgraph'
args['bedLog'] = args['bedPrf'] + '.bedprocess.log'


# Generate object to store commands
moabJobs = moab.moabJobs()
# Create trim command and add to job dictionary
trimCommand = fastqTrim.cutadaptTrimPaired(read1In = args['read1'],
    read2In = args['read2'], read1Out = args['trmRd1'],
    read2Out = args['trmRd2'], path = paths['cutadapt'], length = 25,
    quality = 20, adapter = 'CTGTCTCTTATACACATCT')
trimCommandID = moabJobs.add(trimCommand, stdout = args['trmLog'],
    stderr = args['trmLog'])
# Generate QC commands
read1Fastqc = fastqQC.fastQC(inFile = args['trmRd1'], outDir = args['smpFsq'],
        path = paths['fastqc'])
read2Fastqc = fastqQC.fastQC(inFile = args['trmRd2'], outDir = args['smpFsq'],
        path = paths['fastqc'])
# Combine FASTQ QC command and add to job dictionary
fastqcCommand = '%s && %s' %(read1Fastqc, read2Fastqc)
fastqcCommandID = moabJobs.add(fastqcCommand, stdout = args['fqcLog'],
    stderr = args['fqcLog'], dependency = [trimCommandID])
# Generate command to perform alignment and add to job dictionary
alignCommand = fastqAlign.bowtie2Align(index = args['<bwt2index>'],
    outFile = args['alnBam'], read1 = args['read1'], read2 = args['read2'],
    bowtie2Path = paths['bowtie2'], threads = args['--threads'],
    readGroup = args['name'], sampleName = args['name'],
    libraryID = args['prefix'], discordant = False, mixed = False,
    maxInsert = 2000, check = True, samtoolsPath = paths['samtools'])
alignCommandID = moabJobs.add(alignCommand, stdout = args['alnLog'],
    stderr = args['alnLog'], processors = args['--threads'],
    dependency = [trimCommandID])
# Generate command to deduplicate BAM and add to job dictionary
dedupCommand = picard.markDuplicates(inBam = args['alnBam'],
    outBam = args['mdpBam'], logFile = args['mdpLog1'],
    picardPath = paths['picard'], javaPath = paths['java'],
    removeDuplicates = False, delete = True)
dedupCommandID = moabJobs.add(dedupCommand, stdout = args['mdpLog2'],
    stderr = args['mdpLog2'], dependency = [alignCommandID])
# Generate command to create insert size metrics
insertCommand = picard.collectInsertSizeMetrics(inBam = args['mdpBam'],
    outPrefix = args['insPrf'], picardPath = paths['picard'],
    javaPath = paths['java'])
insertCommandID = moabJobs.add(insertCommand, stdout = args['insLog'],
    stderr = args['insLog'], dependency = [dedupCommandID])
# Generate command to create bed file and add to job dictionary
bam2bedCommand = '%s %s %s %s --minMapQ %s --size %s' %(
    paths['python'], paths['bam2bed'], args['mdpBam'], args['intBed'],
    args['--minMapQ'], args['--interval'])
if args['--rmDup']:
    bam2bedCommand += ' --rmDup'
bam2bedCommandID = moabJobs.add(bam2bedCommand, stdout = args['intLog'],
    stderr = args['intLog'], dependency = [dedupCommandID])
# Generate command to sort bed file and generate bedgraph
sortCommand = bedtools.sortBed(inBed = args['intBed'], outBed = args['srtBed'],
    path = paths['bedtools'], delete = True)
bgCommand = bedtools.bed2bedGraph(inBed = args['srtBed'], outBG = args['bedGrp'],
    chrFile = args['<chrfile>'], path = paths['bedtools'], delete = False)
# Join sort and bedgraph commands and add to job dictionary
jointCommand = '%s && %s' %(sortCommand, bgCommand)
joinCommandID = moabJobs.add(jointCommand, stdout = args['bedLog'],
    stderr = args['bedLog'], dependency = [bam2bedCommandID])
# Submit jobs
moabJobs.submit(verbose = True)
