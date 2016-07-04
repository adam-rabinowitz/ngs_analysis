'''annotateVariants.py
    
    Usage:
    
    varscan.py <sampledata> <bamdir> <varscandir> <outdir> <annovar>
        <buildver> <database> <fasta>
    
'''
# Import modules
import os
from general_python import moab, docopt, toolbox
from ngs_python.variant import varscan
from ngs_python.bam import pysamfunc
# Extract arguments
args = docopt.docopt(__doc__, version = 'v1')

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
# Check each idividual has a normal
for individual in sampleDict:
    # Check for normal
    if not 'NORM' in sampleDict[individual]:
        raise IOError('no normal for %s' %(individual))

# Extract BAM files from input directory
bamFiles = []
for obj in os.listdir(args['<bamdir>']):
    obj = os.path.join(args['<bamdir>'], obj)
    if os.path.isfile(obj) and obj.endswith('.bam'):
        bamFiles.append(obj)
# Associate BAM files with individual samples
for individual in sampleDict:
    for tissue in sampleDict[individual]:
        name = sampleDict[individual][tissue][0]
        prefix = os.path.join(args['<bamdir>'], name)
        for bam in bamFiles:
            if bam.startswith(prefix):
                sampleDict[individual][tissue].append(bam)
# Check each idividual has a normal and every tissue has a BAM
for individual in sampleDict:
    # check for individual BAM files
    for tissue in sampleDict[individual]:
        if len(sampleDict[individual][tissue]) < 2:
            raise IOError('no bam found for %s %s' %(individual, tissue))
        if len(sampleDict[individual][tissue]) > 2:
            raise IOError('multiple bams found for %s %s' %(individual,
                tissue))

# Extract snp files from varscan directory
snpFiles = []
for obj in os.listdir(args['<varscandir>']):
    obj = os.path.join(args['<varscandir>'], obj)
    if os.path.isfile(obj) and obj.endswith('.snp'):
        snpFiles.append(obj)
# Associate snp files with individual samples
for individual in sampleDict:
    for tissue in sampleDict[individual]:
        name = sampleDict[individual][tissue][0]
        prefix = os.path.join(args['<varscandir>'], name)
        for snp in snpFiles:
            if snp.startswith(prefix):
                sampleDict[individual][tissue].append(snp)
# Check each idividual has a single indel file
for individual in sampleDict:
    # check for individual  files
    for tissue in sampleDict[individual]:
        if tissue == 'NORM':
            continue
        if len(sampleDict[individual][tissue]) < 3:
            raise IOError('no snp file found for %s %s' %(individual, tissue))
        if len(sampleDict[individual][tissue]) > 3:
            raise IOError('multiple snp files found for %s %s' %(individual,
                tissue))

# Extract indel files from varscan directory
indelFiles = []
for obj in os.listdir(args['<varscandir>']):
    obj = os.path.join(args['<varscandir>'], obj)
    if os.path.isfile(obj) and obj.endswith('.indel'):
        indelFiles.append(obj)
# Associate indel files with individual samples
for individual in sampleDict:
    for tissue in sampleDict[individual]:
        name = sampleDict[individual][tissue][0]
        prefix = os.path.join(args['<varscandir>'], name)
        for indel in indelFiles:
            if indel.startswith(prefix):
                sampleDict[individual][tissue].append(indel)
# Check each idividual has a single indel file
for individual in sampleDict:
    # check for individual BAM files
    for tissue in sampleDict[individual]:
        if tissue == 'NORM':
            continue
        if len(sampleDict[individual][tissue]) < 4:
            raise IOError('no indel file found for %s %s' %(individual,
                tissue))
        if len(sampleDict[individual][tissue]) > 4:
            raise IOError('multiple indel files found for %s %s' %(individual,
                tissue))

for individual in sampleDict:
    # Create variables to stroe data
    sampleNames = ['NORM']
    bamFiles = []
    variants = set()
    print '\n{}'.format(individual)
    # Extract bam for normal
    normName, normBam = sampleDict[individual].pop('NORM')
    bamFiles.append(normBam)
    # Create and store varsan commands
    for tissue, tissueData in sorted(sampleDict[individual].items()):
        # Store sample names and bam files
        sampleNames.append(tissue)
        bamFiles.append(tissueData[1])
        # Extract snp and indel variants
        snpSet = varscan.variantSet(tissueData[2], somatic = True)
        print '{} snp {}'.format(tissue, len(snpSet))
        variants.update(snpSet)
        indelSet = varscan.variantSet(tissueData[3], somatic = True)
        print '{} indel {}'.format(tissue, len(indelSet))
        variants.update(indelSet)
    # Convert variant set to list
    variants = list(variants)
    print 'total somatic {}'.format(len(variants))
    # Create filenames
    outFile = os.path.join(args['<outdir>'], individual + '.annotation')
    tempPrefix = os.path.join(args['<outdir>'], individual + '.tmp')
    varData = pysamfunc.calculateVariantMetrics(
        variantList = variants, bamList = bamFiles, sampleNames = sampleNames,
        annovarPath = args['<annovar>'], buildver = args['<buildver>'],
        database = args['<database>'], tempprefix = tempPrefix,
        minMapQ = 20, minBaseQ = 20, groupdel = True, homo = True,
        complexity = True, altQualNormal = 10, fasta = args['<fasta>']
    )
    varData.to_csv(outFile, sep = '\t', index = False)
