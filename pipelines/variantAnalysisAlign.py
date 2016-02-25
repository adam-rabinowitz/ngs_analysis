'''alignSortIndex.py

Usage:
    
    alignSortIndex.py <sampledata> <indir> <outdir> <index> <indelvcf>
        <pathfile> [--threads=<threads>]
    
Options:
    
    --threads=<threads>  Number of threads [default: 1]
    --pair=<pair>        Input reads should be paired
    --help               Output this message
    
'''
# Import required modules
import os
from general_python import docopt
from ngs_python.fastq import fastqFind, fastqAlign
from ngs_python.bam import picard, gatk
from general_python import moab
# Extract arguments
args = docopt.docopt(__doc__, version = 'v1')
# Check numerical arguments
args['--threads'] = int(args['--threads'])
# Split sample data
args['prefix'], args['name'] = args['<sampledata>'].split(',')
# Split input directories into a list
args['<indir>'] = args['<indir>'].split(',')
# Read in path file
paths = {}
with open(args['<pathfile>'], 'r') as pfile:
    for line in pfile:
        program, path = line.strip().split('\t')
        paths[program] = path
print paths
# Create folder for log files
args['logdir'] = os.path.join(args['<outdir>'], 'logData')
if not os.path.isdir(args['logdir']):
    os.mkdir(args['logdir'])
# Find fastq files
read1, read2 = fastqFind.findFastq(
    prefix = args['prefix'], dirList = args['<indir>'], pair = True,
    gzip = True
)
if len(read1) != 1 and len(read2) != 1:
    raise IOError('Failure to find single paired FASTQ files')
# Generate output files
bamPrefix = os.path.join(args['<outdir>'], args['name'])
logPrefix = os.path.join(args['logdir'], args['name'])
outfiles = {
    'initialbam' : bamPrefix + '.bam',
    'dedupbam' : bamPrefix + '_dedup.bam',
    'realignbam' : bamPrefix + '_dedup_realign.bam',
    'listfile' : logPrefix + '_target.list',
    'alignlog' : logPrefix + '_align.log',
    'deduplog1' : logPrefix + '_dedup_1.log',
    'deduplog2' : logPrefix + '_dedup_2.log',
    'realignlog' : logPrefix + '_realign.log'
}
# Generate command for alignment
alignCommand = fastqAlign.bwaMemAlign(
    index = args['<index>'], outFile = outfiles['initialbam'],
    read1 = read1[0], read2 = read2[0], bwaPath = paths['bwa'],
    threads = args['--threads'], sampleName = args['name'],
    libraryID = args['prefix'], readGroup = 1, markSecondary = True,
    checkIndex = True, samtoolsPath = paths['samtools'], memory = 2,
    nameSort = False
)
# Mark duplicates using picard
dedupCommand = picard.markDuplicates(
    inBam = outfiles['initialbam'], outBam = outfiles['dedupbam'],
    logFile = outfiles['deduplog1'], picardPath = paths['picard'],
    javaPath = paths['java'], removeDuplicates = False, delete = True
)
# Perform local realignment
realignCommand = gatk.gatkRealign(
    inBam = outfiles['dedupbam'], outBam = outfiles['realignbam'],
    inVcf = args['<indelvcf>'], reference = args['<index>'],
    javaPath = paths['java'], gatkPath = paths['gatk'], delete = True,
    threads = 4, listFile = outfiles['listfile']
)
# Submit jobs
alignJobID = moab.submitjob(
    alignCommand, stdout = outfiles['alignlog'], stderr = outfiles['alignlog'],
    processor = args['--threads']
)
print alignCommand 
print alignJobID
dedupJobID = moab.submitjob(
    dedupCommand, stdout = outfiles['deduplog2'],
    stderr = outfiles['deduplog2'], processor = 1,
    dependency = alignJobID
)
print dedupCommand
print dedupJobID
realignJobID = moab.submitjob(
    realignCommand, stdout = outfiles['realignlog'],
    stderr = outfiles['realignlog'], processor = args['--threads'],
    dependency = dedupJobID
)
print realignCommand
print realignJobID
