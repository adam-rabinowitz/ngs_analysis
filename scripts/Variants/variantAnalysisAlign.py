'''variantAnalysisAlign.py

Usage:
    
    variantAnalysisAlign.py <sampledata> <indir> <outdir> <index>
        <indelvcf> <snpvcf> <pathfile> [--threads=<threads>]
    
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
paths ={}
with open(args['<pathfile>'], 'r') as pfile:
    for line in pfile:
        program, path = line.strip().split('\t')
        paths[program] = path
# Create folder for log files
args['logdir'] = os.path.join(args['<outdir>'], args['name'] + '_log')
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
    'recalbam' : bamPrefix + '_dedup_realign_recal.bam',
    'listfile' : logPrefix + '_target.list',
    'bsqrfile' : logPrefix + '_bsqr.grp',
    'alignlog' : logPrefix + '_align.log',
    'deduplog1' : logPrefix + '_dedup_1.log',
    'deduplog2' : logPrefix + '_dedup_2.log',
    'realignlog' : logPrefix + '_realign.log',
    'recallog' : logPrefix + '_recal.log'
}
# Generate command for alignment
alignCommand = fastqAlign.bwaMemAlign(
    index = args['<index>'], outFile = outfiles['initialbam'],
    read1 = read1[0], read2 = read2[0], bwaPath = paths['bwa'],
    threads = args['--threads'], sampleName = args['name'],
    libraryID = args['prefix'], readGroup = 1, platform = 'ILLUMINA',
    markSecondary = True, check = True, samtoolsPath = paths['samtools'],
    memory = 2, nameSort = False
)
# Mark duplicates using picard
dedupCommand = picard.markDuplicates(
    inBam = outfiles['initialbam'], outBam = outfiles['dedupbam'],
    logFile = outfiles['deduplog1'], picardPath = paths['picard'],
    javaPath = paths['java'], removeDuplicates = True, delete = True
)
# Perform local realignment
realignCommand = gatk.gatkRealign(
    inBam = outfiles['dedupbam'], outBam = outfiles['realignbam'],
    inVcf = args['<indelvcf>'], reference = args['<index>'],
    javaPath = paths['java'], gatkPath = paths['gatk'], delete = True,
    threads = 4, listFile = outfiles['listfile']
)
recalCommand = gatk.bsqr(
    inBam = outfiles['realignbam'], outBam = outfiles['recalbam'],
    inVcf = args['<snpvcf>'], reference = args['<index>'],
    bsqrTable = outfiles['bsqrfile'], javaPath = paths['java'],
    gatkPath = paths['gatk'], delete = True
)
# Add jobs to moabJobs object
moabJobs = moab.moabJobs()
alignID = moabJobs.add(
    alignCommand, stdout = outfiles['alignlog'], stderr = outfiles['alignlog'],
    processors = args['--threads']
)
dedupID = moabJobs.add(
    dedupCommand, stdout = outfiles['deduplog2'],
    stderr = outfiles['deduplog2'], processors = 1,
    dependency = [alignID]
)
realignID = moabJobs.add(
    realignCommand, stdout = outfiles['realignlog'],
    stderr = outfiles['realignlog'], processors = args['--threads'],
    dependency = [dedupID]
)
recalID = moabJobs.add(
    recalCommand, stdout = outfiles['recallog'],
    stderr = outfiles['recallog'], processors = 1,
    dependency = [realignID]
)
# Submit jobs and print moab identifiers
moabJobs.submit(verbose = True)
