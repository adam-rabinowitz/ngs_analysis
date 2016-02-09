################################################################################
"""VariantAnalysis.py

Usage:
    
    VariantAnalysis.py <sampledata> <indir> <outdir> <genprefix> <indelvcf>
        [--baseq=<trimq>] [--alignq=<alignq>] [--threads=<threads>]
    
Options:
    
    --baseq=<label>      Minimum base quality [default: 20]
    --alignq=<alignq>    Minimum alignment quality [default: 20]
    --threads=<threads>  Number of threads [default: 4]
    --help               Output this message
    
"""
## Initialise analysis
################################################################################
# Import functions from external modules
import subprocess
import sys
import os
import collections
import operator
# Import custom functions
from general_functions import docopt
from ngs_python.fastq import fastqFind, fastqAlign
from ngs_python.bam import samtools
# Extract arguments and check numerical values
args = docopt.docopt(__doc__,version = 'v1')
args['--baseq'] = int(args['--baseq'])
args['--alignq'] = int(args['--alignq'])
args['--threads'] = int(args['--threads'])
# Extract sample name and prefix
args['sample'], args['prefix'] = args['<sampledata>'].split(',')
# Extract input fastq files
args['read1'], args['read2'] = fastqFind.findIlluminaFastq(
    prefix = args['prefix'],
    dirList = args['<indir>'].split(','),
    pair = True
)
# Generate output directories
args.alignDir = args['<outdir>'] + "bamFiles/"
args.varscanDir = args['<outdir>'] + "varscan/"
for folder in [args.alignDir, args.varscanDir]:
    if not os.path.exists(folder):
        os.mkdir(folder)
# Set paths
paths = {
    'bwa' : '/farm/babs/redhat6/software/bwa-0.7.12/bwa',
    'samtools' : '/farm/babs/redhat6/software/samtools-1.2/bin/samtools'
}

################################################################################
## Perform alignmnets
################################################################################
# Create file names
outPrefix = args['<outdir>'] + args['sample']
files = {
    'initialSam' : outPrefix + '_initial.sam',
    'initialBam' : outPrefix + '_initial.bam',
    'sortBam' : outPrefix + '_sort.bam',
    'dedupBam' : outPrefix + '_dedup.bam',
    'dedupLog' : outPrefix + '_dedup.log',
    'targetIntervals' : outPrefix + '_target.interval_list',
    'realignBam' : outPrefix + '_dedup_realign.bam',
    'pileup' : outPrefix + '.pileup',
    'pileup_nz' : outPrefix + '_nonzero.pileup',
    'processLog' : outPrefix + '.log'
}
alignSortCommand = []
# Extract alignment command
alignSortCommand.append(
    fastqAlign.bwaMemAlign(
        index = args['<genprefix>'],
        outSam = files['initialSam'],
        read1 = args['read1'][0],
        read2 = args['read2'][0],
        threads = args['--threads'],
        markSecondary = False,
        path = paths['bwa']
    )
)
# Convert SAM to BAM
alignSortCommand.append(
    samtools.sort(
        inFile = files['initialSam'],
        outFile = files['sortBam'],
        threads = args['--threads'],
        delete = True,
        path = paths['samtools']
    )
)
# Join commands
alignSortCommand = ' && '.join(alignSortCommand)
print alignSortCommand
#		inSam = initialSam,
#		outBam = initialBam,
#		execute = False
#	)
#	# Sort BAM file
#	sortBamCommand = alignment_functions.samtools_sort_bam(
#		inBam = initialBam,
#		outBam = sortBam,
#		threads = 4,
#		execute = False
#	)
#	# Index bam
#	# Mark duplicates in bam
#	markDuplicateCommand = alignment_functions.picard_mark_duplicates(
#		inBam = sortBam,
#		outBam = dedupBam,
#		logFile = dedupLog,
#		remove_duplicates = True,
#		execute = False
#)
#	# Index bam
#	indexDedupCommand = alignment_functions.samtools_index(
#		inBam = dedupBam,
#		execute = False
#	)
#	# Find intervals for realignment
#	findIntervalsCommand = alignment_functions.gatk_realign_target_creator(
#		inBam = dedupBam,
#		inVcf = indelVcf,
#		reference = args.indexedReference,
#		outList = targetIntervals,
#		threads = 4,
#		execute = False
#	)
#	# Peform local realignment
#	localRealignCommand = alignment_functions.gatk_realign(
#		inBam = dedupBam,
#		inVcf = args.indelVcf,
#		reference = args.indexedReference,
#		target_intervals = targetIntervals,
#		outBam = realignBam,
#		execute = False
#	)
#	# Index bam
#	indexRealignCommand = alignment_functions.samtools_index(
#		inBam = dedupBam,
#		execute = False
#	)
#	# Perform pileup
#	pileupCommand = ['/farm/babs/redhat6/software/samtools-0.1.19/samtools',
#		'mpileup', '-A', '-B', '-d', '10000', '-q', '20', '-Q', '20', '-f',
#		args.indexedReference, realignBam, '>', pileup]
#	# Remove zero pileups
#	removeZeroCommand = ['awk', "'$4 > 0'", pileup, '>', pileup_nz]
#	# Merge commands
#	combinedCommand = '%s > %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; rm %s' %(
#		' '.join(bwaAlignCommand),
#		initialSam,
#		' '.join(samBamCommand),
#		' '.join(sortBamCommand),
#		' '.join(indexInitialCommand),
#		' '.join(markDuplicateCommand),
#		' '.join(indexDedupCommand),
#		' '.join(findIntervalsCommand),
#		' '.join(localRealignCommand),
#		' '.join(indexRealignCommand),
#		' '.join(pileupCommand),
#		' '.join(removeZeroCommand),
#		' '.join([
#				initialSam,
#				initialBam,
#				sortBam + '*',
#				dedupBam + '*',
#				dedupLog,
#				targetIntervals,
#				pileup
#			])
#	)
#	# Submit commands to farm
#	jobID = system_functions.submit_moab_job(
#		command = combinedCommand,
#		stdout = processLog,
#		stderr = processLog,
#		processor = 4
#	)
#	return(jobID)
#
#################################################################################
### Process sample data
#################################################################################
## Create dictionary
#sampleDict = collections.OrderedDict()
## Read through input file and add to dictionary
#sampleFile = open(args.sampleData,'r')
#header = sampleFile.readline().strip().split('\t')
#for line in sampleFile:
#	# Find animal for sample
#	lineData = line.strip().split('\t')
#	animal = lineData[header.index('animal')]
#	if animal in sampleDict:
#		sampleDict[animal].append(lineData)
#	else:
#		sampleDict[animal] = [lineData]
#sampleFile.close()
## Check that each animal has one norm and reorder animal list
#for animal in sampleDict:
#	# Extract tissues
#	tissues = map(
#		operator.itemgetter(header.index('tissue')),
#		sampleDict[animal]
#	)
#	# Check a single NORM sample is present
#	if tissues.count('NORM') != 1:
#		sys.exit('Animal %s does not have a single NORM tissue' %(animal))
#	# Reorder animal list so NORM is first
#	normIndex = tissues.index('NORM')
#	sampleDict[animal] = [sampleDict[animal].pop(normIndex)] + sampleDict[animal]
#
#################################################################################
### Generate alignments
#################################################################################
## Create dictionary to store alignnment job ID
#alignJobDict = collections.OrderedDict()
## Generate alignments for each animal
#for animal in sampleDict:
#	animalJobID = []
#	# Loop through each sample from each animal
#	for sample in sampleDict[animal]:
#		# Extract data for sample
#		ID = sample[header.index('id')]
#		name = sample[header.index('name')]
#		group = sample[header.index('read_group')]
#		tissue = sample[header.index('tissue')]
#		directory = sample[header.index('directory')]
#		reads = system_functions.find_file_prefix(
#			directory,
#			ID
#		)
#		# Check that two reads are present for each sample
#		if len(reads) != 2:
#			sys.exit('Two reads were not found for sample %s' %(name))
#		# Create command to peform alignment
#		jobID = performAlignment(
#			read1 = directory + reads[0],
#			read2 = directory + reads[1],
#			outDir = args.alignDir,
#			sampleName = name,
#			libraryID = ID,
#			readGroup = group,
#			indelVcf = args.indelVcf
#		)
#		# Append job ID
#		animalJobID.append(jobID)
#	# Store job ids for animal
#	alignJobDict[animal] = animalJobID
#
#################################################################################
### Perform varscan copynumber analysis
#################################################################################
## Loop through animals
#for animal in sampleDict:
#	# Extract pileupfile for norm
#	normData = sampleDict[animal][0]
#	normName = normData[header.index('name')]
#	normPileup = args.alignDir + normName + '_nonzero.pileup'
#	# Process tumour samples
#	for sample in sampleDict[animal][1:]:
#		# Extract name of intput file
#		tumourName = sample[header.index('name')]
#		tumourPileup = args.alignDir + tumourName + '_nonzero.pileup'
#		# Create somatic command
#		somaticCommand = variant_functions.varscan_somatic(
#			normal_pileup = normPileup,
#			tumour_pileup = tumourPileup,
#			outPrefix = args.varscanDir + tumourName,
#			execute = False
#		)
#		# Create copynumber command
#		copyCommand = variant_functions.varscan_copynumber(
#			normal_pileup = normPileup,
#			tumour_pileup = tumourPileup,
#			outPrefix = args.varscanDir + tumourName,
#			execute = False
#		)
#		# Create copycaller command
#		callerCommand = variant_functions.varscan_copycaller(
#			copynumber = args.varscanDir + tumourName + '.copynumber',
#			output = args.varscanDir + tumourName + '.copycalled',
#			execute = False
#		)
#		# Merge commands and submit
#		combinedCommand = '%s' %(
#			' '.join(somaticCommand)
#		)
#		# Submit varscan copy with dependencies
#		jobID = system_functions.submit_moab_job(
#			combinedCommand,
#			stdout = args.varscanDir + tumourName + '.log',
#			stderr = args.varscanDir + tumourName + '.log',
#			processor = 1,
#			dependency = alignJobDict[animal]
#		)
#		print '%s: %s' %(tumourName, jobID)
