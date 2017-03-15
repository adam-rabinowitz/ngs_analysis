'''bowtie2Align.py

Usage:
    
    bowtie2Align.py <prefix> <indir> <index> <outbam> <modules>
        [--unpaired] [--namesort] [--samplename=<sn>] [--libraryid=<li>]
        [--readgroup=<rg>] [--platform=<pl>] [--trim=<tm>] [--maxinsert=<mi>]
    
Options:
    
    --unpaired         Input reads are unpaired.
    --namesort         Name sort output file.
    --samplename=<sm>  Sample name [default: Unknown].
    --libraryid=<li>   Library id [default: Unknown].
    --platform=<pl>    Platform [default: Unknown].
    --readgroup=<rg>   Read group [default: 1].
    --trim=<tm>        Sequence to trim from reads.
    --maxinsert=<mi>   Maximum inset size.
    --help             Output this message.
    
'''
# Import required modules
from general_python import docopt
from ngs_python.fastq import fastqFind, fastqAlign, fastqTrim, fastqQC
from ngs_python.bam import picard
from general_python import slurm
import itertools
import os
import subprocess
# Extract and check arguments
args = docopt.docopt(__doc__, version = 'v1')
args['<prefix>'] = args['<prefix>'].split(',')
args['<indir>'] = args['<indir>'].split(',')
args['<indir>'] = [os.path.abspath(x) for x in args['<indir>']]
if args['--maxinsert'] is not None:
    args['--maxinsert'] = int(args['--maxinsert'])
if not args['<outbam>'].endswith('.bam'):
    raise ValueError("Output file must end with '.bam'")
args['<outbam>'] = os.path.abspath(args['<outbam>'])
args['<outprefix>'] = args['<logfolder>'] = args['<outbam>'][:-4]
args['<logprefix>'] = os.path.join(
    args['<logfolder>'], os.path.basename(args['<logfolder>'])
)
if not os.path.isdir(args['<logfolder>']):
    os.mkdir(args['<logfolder>'])
# Parse module file and create job dictionary
pmDict = slurm.parsePathModule(args['<modules>'])
slurmJobs = slurm.submitJobs()

# Find fastq files
read1List = []
read2List = []
if args['--unpaired']:
    for prefix in args['<prefix>']:
        read1 = fastqFind.findFastq(
            prefix=prefix, dirList=args['<indir>'], pair=False)
        read1List.extend(read1)
    read1List = [os.path.abspath(x) for x in read1List]
else:
    for prefix in args['<prefix>']:
        read1, read2 = fastqFind.findFastq(
            prefix=prefix, dirList=args['<indir>'], pair=True)
        read1List.extend(read1)
        read2List.extend(read2)
# Set absolute path to fastq files and check names
read1List = [os.path.abspath(x) for x in read1List]
read2List = [os.path.abspath(x) for x in read2List]
if len(read1List) == 0:
    raise IOError('Failed to find FASTQ files')
for f in read1List + read2List:
    if not f.endswith('.fastq.gz'):
        raise IOError('Unrecognised FASTQ filename format: {}'.format(f))

# Create trim commands if required
trimJobList = []
qcJobList = []
if args['--trim']:
    # Create output file lists and loop through input files
    trim1List = []
    trim2List = []
    for count, (read1, read2) in enumerate(
            itertools.izip_longest(read1List, read2List, fillvalue=None)
        ):
        # Create and store output FASTQ files
        count += 1
        trim1 = args['<outprefix>'] + '.trim{}.R1.fastq.gz'.format(count)
        trim1List.append(trim1)
        if read2:
            trim2 = args['<outprefix>'] + '.trim{}.R2.fastq.gz'.format(count)
            trim2List.append(trim2)
        else:
            trim2 = None
        # Create output log file
        logFile = args['<logprefix>'] + '.trim{}.log'.format(count)
        # Create and store trim command
        trimCommand = fastqTrim.cutadapt(
            read1In=read1, read1Out=trim1, read2In=read2, read2Out=trim2,
            quality=0, adapter=args['--trim'], length = 25,
            path=pmDict[('cutadapt', 'path')]
        )
        trimJobID = slurmJobs.add(
            trimCommand, processors=1, memory=6, stdout=logFile,
            stderr=logFile, modules=pmDict[('cutadapt', 'modules')]
        )
        trimJobList.append(trimJobID)
        # Create fastqc command
        qcLog = args['<logprefix>'] + '.fastqc{}.log'.format(count)
        qcCommand = fastqQC.fastQC(
            inFile=trim1, outDir=args['<logfolder>'],
            path=pmDict[('fastqc', 'path')]
        )
        if trim2:
            qcCommand += ' && {}'.format(
                fastqQC.fastQC(
                    inFile=trim2, outDir=args['<logfolder>'],
                    path=pmDict[('fastqc', 'path')]
                )
            )
        qcJobID = slurmJobs.add(
            qcCommand, processors=1, memory=6, stdout=qcLog, stderr=qcLog,
            modules=pmDict[('fastqc', 'modules')], depend=[trimJobID]
        )
        qcJobList.append(qcJobID)
    # Switch input file names
    read1List = trim1List
    read2List = trim2List

# Perform alignment for reads
unmarkedBAM = args['<outprefix>'] + '.temp.bam'
bowtie2Log = args['<logprefix>'] + '.bowtie2.log'
bowtie2Command = fastqAlign.bowtie2Align(
    index=args['<index>'], outFile=unmarkedBAM, read1=read1List,
    read2=read2List, bowtie2Path=pmDict[('bowtie2', 'path')],
    threads=16, readGroup=args['--readgroup'],
    sampleName=args['--samplename'], libraryID=args['--libraryid'],
    platform=args['--platform'], discordant=False, mixed=False,
    maxInsert=args['--maxinsert'], check=True,
    samtoolsPath=pmDict[('samtools', 'path')],
    memory=5, nameSort = False
)
bowtie2JobID = slurmJobs.add(
    bowtie2Command, processors=16, memory=6, stdout=bowtie2Log,
    modules=pmDict[('bowtie2', 'modules')] + pmDict[('samtools', 'modules')],
    depend=trimJobList, stderr=bowtie2Log
)

# Delete trim files if generated
if args['--trim']:
    rmCommand = 'rm {}'.format(' '.join(read1List + read2List))
    slurmJobs.add(rmCommand, depend=qcJobList + [bowtie2JobID])

# Mark duplicates in BAM files
dedupLog = args['<logprefix>'] + '.dedup.log'
dedupCommand = picard.markDuplicates(
    inBam=unmarkedBAM, outBam=args['<outbam>'], logFile=dedupLog + '1',
    picardPath=pmDict[('picard', 'path')], javaPath='java',
    removeDuplicates=False, delete=True, memory=30
)
# Add job to queue
dedupJobID = slurmJobs.add(
    command=dedupCommand, processors=6, memory=6, stdout=dedupLog + '2',
    stderr=dedupLog + '2', modules=pmDict[('picard', 'modules')],
    depend=[bowtie2JobID]
)

#Submit jobs
slurmJobs.submit(verbose=True, check_sub=False)
