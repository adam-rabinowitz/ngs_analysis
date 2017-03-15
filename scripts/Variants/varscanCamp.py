'''varscan.py
    
    Usage:
    
    varscan.py <sampledata> <indir> <outdir> <reference> <paths>
    
'''
# Import modules
import os
from general_python import moab, docopt, toolbox
from ngs_python.variant import varscan
from ngs_python.bam import samtools
from general_python import slurm
# Extract arguments
args = docopt.docopt(__doc__, version = 'v1')
args['<indir>'] = os.path.abspath(args['<indir>'])
args['<outdir>'] = os.path.abspath(args['<outdir>'])
args['<reference>'] = os.path.abspath(args['<reference>'])
args['<paths>'] = os.path.abspath(args['<paths>'])
# Extract paths and modules
pmdata = slurm.parsePathModule(args['<paths>'])
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
slurmJobs = slurm.submitJobs()
for individual, indvData in sampleDict.items():
    # Create lists to store job IDs
    pileupFileList = []
    pileupJobList = []
    varscanJobList = []
    # Create log directory
    logDir = os.path.join(args['<outdir>'], individual + '.log')
    if not os.path.isdir(logDir):
        os.mkdir(logDir)
    # Create and store pileup commands
    for tissue, tissueData in indvData.items():
        # Extract sample name and BAM file
        name, bam = tissueData
        # Create and store file names
        pileup = os.path.join(args['<outdir>'], name + '.mpileup')
        pileupFileList.append(pileup)
        pileupLog = os.path.join(logDir, name + '.mpileup.log')
        # Create pileup command
        pileupCommand = samtools.mpileup(
            inBam = bam, outFile = pileup, reference = args['<reference>'],
            minMapQ = 20, minBaseQ = 20, countOrphans = False,
            countDup = False, disableBAQ = True, maxDepth = 10000,
            minFilter = 1, path = pmdata[('samtools', 'path')]
        )
        # Add pileup command and store jobs
        pileupJobID = slurmJobs.add(
            command = pileupCommand, stdout = pileupLog, stderr = pileupLog,
            memory = 6, modules = pmdata[('samtools', 'modules')],
            processors = 1,
        )
        pileupJobList.append(pileupJobID)
    # Extract data for normal
    normName, normBam = indvData.pop('NORM')
    normPileup = os.path.join(args['<outdir>'], normName + '.mpileup')
    # Create and store varsan commands
    for tissue, tissueData in indvData.items():
        # Extract tissue data
        name = tissueData[0]
        # Create file names
        prefix = os.path.join(args['<outdir>'], name)
        tumourPileup = prefix + '.mpileup'
        somaticLog = os.path.join(logDir, name + '.somatic.log')
        copynumberLog = os.path.join(logDir, name + '.copynumber.log')
        # Create somatic command
        somaticCommand = varscan.somatic(
            mpileup1 = normPileup, mpileup2 = tumourPileup,
            outPrefix = prefix, purity = 0.5, minCovNormal = 8,
            minCovTumour = 6, minHetFreq = 0.1, minHomFreq = 0.75,
            normalPurity = 1.0, tumourPurity = 0.5, pValueHet = 0.99,
            pValueSomatic = 0.05, strandFilter = False, javaPath = 'java',
            varscanPath = pmdata[('varscan', 'path')], memory = 10
        )
        # Add somatic command
        somaticJobID = slurmJobs.add(
            command = somaticCommand, stdout = somaticLog, stderr = somaticLog,
            modules = pmdata[('varscan', 'modules')], depend = pileupJobList,
            processors = 2, memory = 6
        )
        # Create copynumber command
        copynumberCommand = varscan.copynumber(
            mpileup1 = normPileup, mpileup2 = tumourPileup, outPrefix = prefix,
            minBaseQ = 20, minMapQ = 20, minCov = 20, minSegSize = 10,
            maxSegSize = 100, pValue = 0.01, javaPath = 'java',
            varscanPath = pmdata[('varscan', 'path')], memory = 10
        )
        # Add copynumber command
        copynumberJobID = slurmJobs.add(
            command = copynumberCommand, stdout = copynumberLog,
            stderr = copynumberLog, depend = pileupJobList,
            modules = pmdata[('varscan', 'modules')], processors = 2,
            memory = 6
        )
        # Store varscan command
        varscanJobList.extend([somaticJobID, copynumberJobID])
    # Create command to delete pileup
    slurmJobs.add(
        command = 'rm %s' %(' '.join(pileupFileList)),
        depend = varscanJobList
    )
# Submit jobs
slurmJobs.submit(verbose = True, check_sub=False)
