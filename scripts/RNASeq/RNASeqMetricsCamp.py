'''RNASeqMetrics.py

Usage:
    
    RNASeqMetrics.py <indir> <outprefix> [--singleend]
    
Options:
    
    --singleend  Samples are single-end sequenced
    
'''
# Import functions
import sys
import argparse
import os
import pandas
import re
import collections
from general_python import docopt
# Extract arguments
args = docopt.docopt(__doc__,version = 'v1')
# Extract directories in all top directories 
dirDict = {}
dirDict['rsemDirList'] = os.listdir(
    os.path.join(args['<indir>'] + 'rsemAlign'))
dirDict['rsemDirList'].sort()
dirDict['tophatDirList'] = os.listdir(
    os.path.join(args['<indir>'] + 'tophat2Align'))
dirDict['tophatDirList'].sort()
dirDict['fastqDirList'] = os.listdir(
    os.path.join(args['<indir>'] + 'trimFastq'))
dirDict['fastqDirList'].sort()
# Set expected counts
if args['--singleend']:
    expect = 1
else:
    expect = 2
# Check that all directories contain the same sub directories
for number, directory in enumerate(dirDict):
    if number:
        if subDir != dirDict[directory]:
            raise IOError('Sub-directories are not identical')
    else:
        subDir = dirDict[directory]
# Create output data frames
expectedCounts = pandas.DataFrame()
tpmCounts = pandas.DataFrame()
fpkmCounts = pandas.DataFrame()
columns = ['initial_reads', 'trim_loss_rate', 'aligned_reads',
    'rsem_align_rate', 'rsem_expressed_genes', 'tophat_align_rate',
    'tophat_exonic_rate', 'tophat_intragenic_rate', 'tophat_intergenic_rate',
    'tophat_duplication_rate', 'tophat_rRNA_rate']
qcMetrics = pandas.DataFrame( columns = columns,
    index = dirDict['rsemDirList'] )

################################################################################
## Process RSEM Transcript data
################################################################################
# Loop through directory and extract rsem data
os.chdir(os.path.join(args['<indir>'], 'rsemAlign'))
for directory in dirDict['rsemDirList']:
    # Generate name of gene results file
    geneFileName = '%s/%s.genes.results' %(
        directory,
        directory
    )
    # Extract transcript counts
    if os.path.isfile(geneFileName):
        geneData = pandas.read_csv(
            geneFileName,
            sep='\t',
            index_col = 0
        )
        expectedCounts[directory] = geneData['expected_count']
        tpmCounts[directory] = geneData['TPM']
        fpkmCounts[directory] = geneData['FPKM']
        # Store number of expressed genes
        qcMetrics.loc[directory, 'rsem_expressed_genes'] = sum(
            geneData['expected_count'] > 0)
    else:
        print "No RSEM transcript '.cnt' for sample %s" %(directory)
    # Generate name of rsem count file
    countFileName = '%s/%s.stat/%s.cnt' %(
        directory,
        directory,
        directory
    )
    # Extract RSEM transcript counts
    if os.path.isfile(countFileName):
        rsemAlignFile = open(countFileName, 'r')
        alignData = rsemAlignFile.readline().strip().split()
        rsemAlignFile.close()
        qcMetrics.loc[directory, 'aligned_reads'] = int(alignData[3])
        qcMetrics.loc[directory, 'rsem_align_rate'] = ( float(alignData[1]) /
            int(alignData[3]))
    else:
        print "No RSEM transcript '.genes.results' file for sample %s" %(
            directory)
    
################################################################################
## Process Tophat Data
################################################################################
# Loop through directory and extract tophat and rnaseqc data
os.chdir(os.path.join(args['<indir>'], 'tophat2Align'))
for directory in dirDict['tophatDirList']:
    # Generate name of alignment file
    alignFileName = '%s/align_summary.txt' %(
        directory
    )
    # Extract alignment metrics
    if os.path.isfile(alignFileName):
        with open(alignFileName, 'r') as tophatFile:
            alignData = tophatFile.read()
        inputValues = map(int, re.findall('Input\s+:\s+(\d+)',alignData))
        mappedValues = map(float, re.findall('Mapped\s+:\s+(\d+)',alignData))
        if (len(inputValues) == expect and len(mappedValues) == expect):
            qcMetrics.loc[directory,'tophat_align_rate'] = (sum(mappedValues) /
                sum(inputValues))
        else:
            raise IOError('Failed to extract tophat2 alignment summary data')
    # Or report their absence
    else:
        print "No Tophat2 'align_summary.txt' file for sample %s" %(directory)
    # Generate rnaseqc file name
    seqcFileName = '%s/metrics.tsv' %(
        directory
    )
    # Check existence of file and extract metrics
    if os.path.isfile(seqcFileName):
        # Read in rnaseqc metrics
        seqcData = pandas.read_csv(
            seqcFileName,
            sep = '\t'
        )
        qcMetrics.loc[directory, 'tophat_duplication_rate'] = seqcData.loc[0,
            'Duplication Rate of Mapped']
        qcMetrics.loc[directory, 'tophat_exonic_rate'] = seqcData.loc[0,
            'Exonic Rate']
        qcMetrics.loc[directory, 'tophat_intragenic_rate'] = seqcData.loc[0,
            'Intragenic Rate']
        qcMetrics.loc[directory, 'tophat_intergenic_rate'] = seqcData.loc[0,
            'Intergenic Rate']
        qcMetrics.loc[directory, 'tophat_rRNA_rate'] = seqcData.loc[0,
            'rRNA rate']
    else:
        print "No RNASeqC 'metrics.tsv' file for sample %s" %(directory)

###############################################################################
## Process FASTQ data
###############################################################################
# Loop through directory and extract fastq data
os.chdir(os.path.join(args['<indir>'], 'trimFastq'))
for directory in dirDict['fastqDirList']:
    # Generate file name of trim log file
    trimFileName = '%s/%s_cutadapt.log' %(
        directory,
        directory
    )
    # Open file and extract data
    with open(trimFileName, 'r') as trimFile:
        trimData = trimFile.read()
    # Extract trimmed read stats
    if args['--singleend']:
        processed = re.search('Total reads processed:\\s+(\\d{1,3}(,\\d{3})*)',
            trimData).group(1)
        tooShort = re.search('Reads that were too short:\\s+(\\d{1,3}(,\\d{3})*)',
            trimData).group(1)
    else:
        processed = re.search('Total read pairs processed:\\s+(\\d{1,3}(,\\d{3})*)',
            trimData).group(1)
        tooShort = re.search('Pairs that were too short:\\s+(\\d{1,3}(,\\d{3})*)',
            trimData).group(1)
    # Process trimmed reads stats and add to output
    processed = re.sub(',','',processed)
    tooShort = re.sub(',','',tooShort)
    qcMetrics.loc[directory, 'initial_reads'] = processed
    qcMetrics.loc[directory, 'trim_loss_rate'] = (int(tooShort) /
            float(processed))

# Save output data to file
expectedCounts.to_csv(
    args['<outprefix>'] + '_gene.exp.txt',
    sep = '\t'
)
tpmCounts.to_csv(
    args['<outprefix>'] + '_gene.tpm.txt',
    sep = '\t'
)
fpkmCounts.to_csv(
    args['<outprefix>'] + '_gene.fpkm.txt',
    sep = '\t'
)
qcMetrics.to_csv(
    args['<outprefix>'] + '.qc.txt',
    sep = '\t'
)
