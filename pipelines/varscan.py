'''varscan.py
    
    Usage:
    
    varscan.py <sampledata> <indir> <outdir> <reference> <paths>
    
'''
# Import modules
import os
from general_python import moab, docopt, toolbox
from ngs_python.variant import varscan
from ngs_python.bam import samtools
# Extract arguments
args = docopt.docopt(__doc__, version = 'v1')
# Extract paths
paths = toolbox.fileDict(args['<paths>'], sep = '\t')
# Read through input file and create sample dictionary
sampleDict = {}
with open(args['<sampledata>'], 'r') as f:
    header = f.readline().strip().split('\t')
    nIndex = header.index('name')
    iIndex = header.index('individual')
    tIndex = header.index('tissue')
    for line in f:
        # Extract sample data
        lineData = line.strip().split('\t')
        individual = lineData[iIndex]
        tissue = lineData[tIndex]
        name = lineData[nIndex]
        # Add to sample dictionary
        if individual in sampleDict:
            sampleDict[individual][tissue] = [name]
        else:
            sampleDict[individual] = {tissue : [name]}
# Extract BAM files from input directory
bamFiles = []
for obj in os.listdir(args['<indir>']):
    obj = os.path.join(args['<indir>'], obj)
    if os.path.isfile(obj) and obj.endswith('.bam'):
        bamFiles.append(obj)
# Associate BAM files with individual samples
for individual in sampleDict:
    for tissue in sampleDict[individual]:
        name = sampleDict[individual][tissue][0]
        prefix = os.path.join(args['<indir>'], name)
        for bam in bamFiles:
            if bam.startswith(prefix):
                sampleDict[individual][tissue].append(bam)
# Check each idividual has a normal and every tissue has a BAM
for individual in sampleDict:
    # Check for normal
    if not 'NORM' in sampleDict[individual]:
        raise IOError('no normal for %s' %(individual))
    # check for individual BAM files
    for tissue in sampleDict[individual]:
        if len(sampleDict[individual][tissue]) < 2:
            raise IOError('no bam found for %s %s' %(individual, tissue))
        if len(sampleDict[individual][tissue]) > 2:
            raise IOError('multiple bams found for %s %s' %(individual,
                tissue))
# Loop through individuals data
for individual, indvData in sampleDict.items():
    # Create objects to store jobs
    moabJobs = moab.moabJobs()
    # Create lists to store job IDs
    pileupFileList = []
    pileupJobList = []
    varscanJobList = []
    # Create log directory
    logDir = os.path.join(args['<outdir>'], individual + '_Log')
    if not os.path.isdir(logDir):
        os.mkdir(logDir)
    # Create and store pileup commands
    for tissue, tissueData in indvData.items():
        # Extract sample name and BAM file
        name, bam = tissueData
        # Create and store file names
        pileup = os.path.join(args['<outdir>'], name + '.mpileup')
        pileupFileList.append(pileup)
        pileupLog = os.path.join(logDir, name + '_mpileup.log')
        # Create pileup command
        pileupCommand = samtools.mpileup(
            inBam = bam, outFile = pileup, reference = args['<reference>'],
            minMapQ = 20, minBaseQ = 20, countOrphans = False,
            countDup = False, disableBAQ = True, maxDepth = 10000,
            minFilter = 1, path = paths['samtools']
        )
        # Add pileup command and store jobs
        pileupJobID = moabJobs.add(
            command = pileupCommand,
            stdout = pileupLog,
            stderr = pileupLog
        )
        pileupJobList.append(pileupJobID)
    # Extract data for normal
    normName, normBam = indvData.pop('NORM')
    normPileup = os.path.join(args['<outdir>'], normName + '.mpileup')
    # Create and store varsan commands
    for tissue, tissueData in indvData.items():
        print individual, tissue
        # Extract tissue data
        name = tissueData[0]
        # Create file names
        prefix = os.path.join(args['<outdir>'], name)
        tumourPileup = prefix + '.mpileup'
        somaticLog = os.path.join(logDir, name + '_somatic.log')
        copynumberLog = os.path.join(logDir, name + '_copynumber.log')
        # Create somatic command
        somaticCommand = varscan.somatic(
            mpileup1 = normPileup, mpileup2 = tumourPileup,
            outPrefix = prefix, purity = 0.5, minCovNormal = 8,
            minCovTumour = 6, minHetFreq = 0.1, minHomFreq = 0.75,
            normalPurity = 1.0, tumourPurity = 0.5, pValueHet = 0.99,
            pValueSomatic = 0.05, strandFilter = False, javaPath = paths['java'],
            varscanPath = paths['varscan']
        )
        # Add somatic command
        somaticJobID = moabJobs.add(
            command = somaticCommand,
            stdout = somaticLog,
            stderr = somaticLog,
            dependency = pileupJobList
        )
        # Create copynumber command
        copynumberCommand = varscan.copynumber(
            mpileup1 = normPileup, mpileup2 = tumourPileup, outPrefix = prefix,
            minBaseQ = 20, minMapQ = 20, minCov = 20, minSegSize = 10,
            maxSegSize = 100, pValue = 0.01, javaPath = paths['java'],
            varscanPath = paths['varscan']
        )
        # Add copynumber command
        copynumberJobID = moabJobs.add(
            command = copynumberCommand,
            stdout = copynumberLog,
            stderr = copynumberLog,
            dependency = pileupJobList
        )
        # Store varscan command
        varscanJobList.extend([somaticJobID, copynumberJobID])
    # Create command to delete pileup
    moabJobs.add(
        command = 'rm %s' %(' '.join(pileupFileList)),
        dependency = varscanJobList
    )
    # Submit jobs
    jobData = moabJobs.submit(verbose = True)
    # Store commands in output file
    moabFile = os.path.join(logDir, individual + '.moab')
    with open(moabFile, 'w') as out:
        for jobID, command, processors, stdout, stderr, depend in jobData:
            out.write('JOB ID: %s\nCOMMAND: %s\nPROCESSORS: %s\nSTDOUT: %s\n'\
                'STDERR: %s\nDEPEND: %s\n\n' %(jobID, command, processors,
                stdout, stderr, depend))
