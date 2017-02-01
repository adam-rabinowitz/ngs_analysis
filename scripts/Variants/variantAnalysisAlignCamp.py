'''variantAnalysisAlign.py

Usage:
    
    variantAnalysisAlign.py <sampledata> <indir> <outdir> <index>
        <snpvcf> <indelvcf> <paths>
    
Options:
    
    --pair=<pair>        Input reads should be paired
    --help               Output this message
    
'''
# Import required modules
import os
from general_python import docopt
from ngs_python.fastq import fastqFind, fastqAlign
from ngs_python.bam import picard, gatk
from general_python import slurm
# Extract and process arguments
args = docopt.docopt(__doc__, version = 'v1')
args['prefix'], args['name'] = args['<sampledata>'].split(',')
args['<indir>'] = args['<indir>'].split(',')
args['<indir>'] = [os.path.abspath(x) for x in args['<indir>']]
args['<outdir>'] = os.path.abspath(args['<outdir>'])
args['<index>'] = os.path.abspath(args['<index>'])
args['<indelvcf>'] = os.path.abspath(args['<indelvcf>'])
args['<snpvcf>'] = os.path.abspath(args['<snpvcf>'])
pmdata = slurm.parsePathModule(os.path.abspath(args['<paths>']))
# Create folder for log files
args['logdir'] = os.path.join(args['<outdir>'], args['name'] + '_log')
if not os.path.isdir(args['logdir']):
    os.mkdir(args['logdir'])
# Find fastq files
print(args['<indir>'])
read1, read2 = fastqFind.findFastq(
    prefix = args['prefix'], dirList = args['<indir>'], pair = True)
if len(read1) != 1 or len(read2) != 1:
    raise IOError('Failure to find single paired FASTQ files')
# Generate output files
bamPrefix = os.path.join(args['<outdir>'], args['name'])
logPrefix = os.path.join(args['logdir'], args['name'])
outfiles = {
    'initialbam' : bamPrefix + '.initial.bam',
    'dedupbam' : bamPrefix + '.dedup.bam',
    'bsqrbam' : bamPrefix + '.bsqr.bam',
    'realignbam' : bamPrefix + '.realign.bam',
    'recalbam' : bamPrefix + 'dedup_realign_recal.bam',
    'alignlog' : logPrefix + '.align.log',
    'deduplog1' : logPrefix + '.dedup1.log',
    'deduplog2' : logPrefix + '.dedup2.log',
    'listfile' : logPrefix + '.target.list',
    'bsqrfile' : logPrefix + '.bsqr.grp',
    'bsqrlog' : logPrefix + '.bsqr.log',
    'realignlog' : logPrefix + '.realign.log',
}
# Create object to store and submit jobs
slurmObject = slurm.submitJobs()
# Generate command for alignment and store
alignCommand = fastqAlign.bwaMemAlign(
    index=args['<index>'], outFile=outfiles['initialbam'],
    read1=read1[0], read2=read2[0], bwaPath=pmdata[('bwa', 'path')],
    threads=16, memory=6, sampleName=args['name'],
    libraryID=args['prefix'], readGroup=1, platform='ILLUMINA',
    markSecondary=True, check=True, nameSort = False,
    samtoolsPath=pmdata[('samtools','path')],
)
alignJobID = slurmObject.add(
    command=alignCommand, processors=16, memory=7, stdout=outfiles['alignlog'],
    stderr=outfiles['alignlog'], modules=pmdata[('bwa', 'modules')] +
    pmdata[('samtools', 'modules')], depend = []
)
# Mark duplicates using picard
dedupCommand = picard.markDuplicates(
    inBam = outfiles['initialbam'], outBam=outfiles['dedupbam'],
    logFile=outfiles['deduplog1'], removeDuplicates=True, delete=True,
    picardPath=pmdata[('picard', 'path')], javaPath='java',
    memory=24
)
dedupJobID = slurmObject.add(
    command=dedupCommand, processors=4, memory=7, stdout=outfiles['deduplog2'],
    stderr=outfiles['deduplog2'], depend=[alignJobID],
    modules=pmdata[('picard', 'modules')]
)
# Perform base quality recalibration
bsqrCommand = gatk.bsqr(
    inBam=outfiles['dedupbam'], outBam=outfiles['bsqrbam'],
    inVcf=args['<snpvcf>'], reference=args['<index>'],
    bsqrTable=outfiles['bsqrfile'], javaPath='java',
    gatkPath=pmdata[('gatk', 'path')], delete=True, memory=24
)
bsqrJobID = slurmObject.add(
    command=bsqrCommand, processors=4, memory=7, stdout=outfiles['bsqrlog'],
    stderr=outfiles['bsqrlog'], depend=[dedupJobID],
    modules=pmdata[('gatk', 'modules')]
)
# Perform local realignment
realignCommand = gatk.gatkRealign(
    inBam=outfiles['bsqrbam'], outBam=outfiles['realignbam'],
    inVcf = args['<indelvcf>'], reference=args['<index>'],
    javaPath='java', gatkPath=pmdata[('gatk', 'path')],
    delete = True, threads = 4, listFile = outfiles['listfile'],
    memory=24
)
realignJobID = slurmObject.add(command=realignCommand, processors=4, memory=7,
    stdout=outfiles['realignlog'], stderr=outfiles['realignlog'],
    depend=[bsqrJobID], modules=pmdata[('gatk', 'modules')]
)
# Submit jobs
slurmObject.submit(verbose=True, check_sub=False)
