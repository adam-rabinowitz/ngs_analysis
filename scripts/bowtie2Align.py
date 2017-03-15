'''bwaAlign.py

Usage:
    
    bwaAlign.py <prefix> <indir> <index> <threads> <outbam> <modules>
        [--unpaired] [--namesort] [--dedupindv] [--samplename=<sn>]
        [--libraryid=<li>] [--platform=<pl>] [--module=<md>]
    
Options:
    
    --unpaired         Input reads are unpaired.
    --namesort         Name sort output file.
    --samplename=<sm>  Sample name [default: Unknown].
    --libraryid=<li>   Library id [default: Unknown].
    --platform=<pl>    Platform [default: Unknown].
    --dedupindv        Dedup individual FASTQ files.
    --help             Output this message.
    
'''
# Import required modules
from general_python import docopt
from ngs_python.fastq import fastqFind, fastqAlign
from ngs_python.bam import picard
from general_python import slurm
import itertools
import os
import subprocess
# Extract and check arguments
args = docopt.docopt(__doc__, version = 'v1')
args['<threads>'] = int(args['<threads>'])
args['<prefix>'] = args['<prefix>'].split(',')
args['<indir>'] = args['<indir>'].split(',')
if not args['<outbam>'].endswith('.bam'):
    raise ValueError("Output file must end with '.bam'")
args['<outbam>'] = os.path.abspath(args['<outbam>'])
# Parse module file
pmDict = slurm.parsePathModule(args['<modules>'])
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
    read1List = [os.path.abspath(x) for x in read1List]
    read2List = [os.path.abspath(x) for x in read2List]
# Raise Error if no Fastq files identified
if len(read1List) == 0:
    raise IOError('Failed to find FASTQ files')
# Create object to store jobs and create log folder
jobObject = slurm.submitJobs()
bamList = []
jobList = []
args['<logfolder>'] = args['<outbam>'][:-4]
os.mkdir(args['<logfolder>'])
# Perform alignment for each pair of FASTQ files
for readgroup, (read1, read2) in enumerate(
        itertools.izip_longest(read1List, read2List)
    ):
    # Modify and format read group and library id
    readgroup = format(readgroup + 1, '03d')
    # Create file names for alignment
    alignLog = os.path.join(args['<logfolder>'], 'rg{}.align.log'.format(
        readgroup))
    alignBam = args['<outbam>'][:-4] + '.rg{}.sort.bam'.format(
        readgroup)
    # Generate command for alignment
    alignCommand = fastqAlign.bwaMemAlign(
        index=os.path.abspath(args['<index>']), outFile=alignBam,
        read1=read1, read2=read2, bwaPath=pmDict[('bwa', 'path')],
        threads=args['<threads>'], sampleName=args['--samplename'],
        libraryID=args['--samplename'], readGroup=readgroup,
        platform=args['--platform'], markSecondary=True, check=True,
        samtoolsPath=pmDict[('samtools', 'path')], memory=6,
        nameSort=False
    )
    # Add job to queue
    alignJobID = jobObject.add(
        command=alignCommand, processors=args['<threads>'],
        modules=pmDict[('bwa', 'modules')] + pmDict[('samtools', 'modules')],
        memory=7, stdout=alignLog, stderr=alignLog
    )
    # Dedup individual files and store outputs
    if args['--dedupindv']:
        # Create file names
        dedupLog = os.path.join(args['<logfolder>'], 'rg{}.dedup.log'.format(
            readgroup))
        dedupBam = args['<outbam>'][:-4] + '.rg{}.dedup.bam'.format(
            readgroup)
        # Create command for marking duplicates
        dedupCommand = picard.markDuplicates(
            inBam=alignBam, outBam=dedupBam, logFile=dedupLog + '1',
            picardPath=pmDict[('picard', 'path')], javaPath='java',
            removeDuplicates=False, delete=True, memory=20
        )
        # Add command to dictionary
        dedupJobID = jobObject.add(
            command=dedupCommand, processors=4,
            modules=pmDict[('picard', 'modules')], memory=6,
            stdout=dedupLog + '2', stderr=dedupLog + '2', depend=[alignJobID]
        )
        # Store output BAM files and job IDs
        bamList.append(dedupBam)
        jobList.append(dedupJobID)
    # Store output for non-dedupped files
    else:
        bamList.append(alignBam)
        jobList.append(alignJobID)
# Merge BAM files
mergedBam = args['<outbam>'][:-4] + '.merged.bam'
if len(bamList) == 1:
    mergeCommand = 'mv {} {}'.format(bamList[0], mergedBam)
    mergeJobID = jobObject.add(mergeCommand, depend=jobList)
else:
    mergeCommand = samtools.merge(
        inBamList=bamList, outBam=mergedBam, delete=False,
        path=pmDict[('samtools', 'path')]
    )
    mergeJobID = jobObject.add(
        mergeCommand, modules=pmDict[('samtools', 'modules')], depend=jobList,
        memory = 24
    )
# Perform final deduplication
if args['--dedupindv']:
    dedupCommand = 'mv {} {}'.format(mergedBam, args['<outbam>'])
    dedupJobID = jobObject.add(dedupCommand)
else:
    # Create file names
    dedupLog = os.path.join(args['<logfolder>'], 'merged.dedup.log')
    # Create command for marking duplicates
    dedupCommand = picard.markDuplicates(
        inBam=mergedBam, outBam=args['<outbam>'], logFile=dedupLog + '1',
        picardPath=pmDict[('picard', 'path')], javaPath='java',
        removeDuplicates=False, delete=True, memory=22
    )
    # Add job to queue
    dedupJobID = jobObject.add(
        command=dedupCommand, processors=1,
        modules=pmDict[('picard', 'modules')], memory=24,
        stdout=dedupLog + '2', stderr=dedupLog + '2',
         depend=[alignJobID]
    )
#Submit jobs
jobObject.submit(verbose=True, check_sub=False)
