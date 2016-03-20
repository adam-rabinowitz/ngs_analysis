'''extractSNP.py
    
    Usage:
    
    extractSNP.py <sampledata> <bamdir> <snpdir> <outdir> <paths>
    
'''
# Import modules
import os
from general_python import docopt, toolbox
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
# Extract BAM files from bam directory
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
# Extract filtered SNP files from snp directory
snpFiles = []
for obj in os.listdir(args['<snpdir>']):
    obj = os.path.join(args['<snpdir>'], obj)
    if os.path.isfile(obj) and obj.endswith('.filter.somatic.snp'):
        snpFiles.append(obj)
# Associate SNP files with individual samples
for individual in sampleDict:
    for tissue in sampleDict[individual]:
        name = sampleDict[individual][tissue][0]
        prefix = os.path.join(args['<snpdir>'], name)
        for snp in snpFiles:
            if snp.startswith(prefix):
                sampleDict[individual][tissue].append(snp)
# Check each non-normal has a SNP file
for individual in sampleDict:
    # check for individual BAM files
    for tissue in sampleDict[individual]:
        if tissue == 'NORM':
            continue
        elif len(sampleDict[individual][tissue]) < 3:
            raise IOError('no SNP file found for %s %s' %(individual, tissue))
        elif len(sampleDict[individual][tissue]) > 3:
            raise IOError('multiple SNP files found for %s %s' %(individual,
                tissue))
# Loop through sample dict and create variant lists and bam lists
bamSNPDict = {}
for individual in sampleDict:
    snpSet = set()
    bamList = []
    for tissue in sampleDict[individual]:
        # Skip normal tissue
        if tissue == 'NORM':
            normBAM = sampleDict[individual][tissue][1]
            continue
        # Extract SNPs for tumours
        else:
            with open(sampleDict[individual][tissue][2], 'r') as snpFile:
                # Skip header
                snpFile.next()
                # Add SNPs to SNP set
                for line in snpFile:
                    chrom, pos, ref, var = line.split('\t')[:4]
                    snpSet.add((chrom, int(pos), ref, var))
    print snpSet
    break
