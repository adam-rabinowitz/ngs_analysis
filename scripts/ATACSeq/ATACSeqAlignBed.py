'''ATACSeqAlignBed.py

Usage:
    
    ATACSeqAlignBed.py <sampledata> <indir> <outdir> <bwt2index>
        <chrfile> <paths> [--minMapQual=<minMapQual>] [--rmDup=<rmDup>]
        [--threads=<threads>]
    
Options:
    
    --minMapQ=<minMapQ>  Minimum mapping quality for reads [default: 20]
    --rmDup=<rmDup>      Flag to ignore duplicate alignments
    --threads=<th>          Number of threads to use [default: 4]

'''
# Import required modules
import os
from ngs_python.fastq import fastqFind, fastqAlign
from general_python import docopt, moab
# Extract and check arguments
args = docopt.docopt(__doc__,version = 'v1')
# Extract sample prefix and name
args['prefix'], args['name'] = args['<sampledata>'].split(',')
# Extract fastq files and check
args['read1'], args['read2'] = fastqFind.findFastq(
    prefix = args['prefix'],
    dirList = args['<indir>'].split(','),
    pair = True)
if len(args['read1']) != len(args['read2']):
    raise IOError('Unequal number of FASTQ files identified')
if len(args['read1']) < 1:
    raise IOError('Insufficient number of FASTQ files identified')
# Create output directories
if not os.path.isdir(outdir):
    raise IOError('Output directory not found')
args['bamDir'] = os.path.join(args['<outdir>'], 'bamFiles')
args['bedDir'] = os.path.join(args['<outdir>'], 'bedFiles')
for directory in [args['bamDir'], args['bedDir']]:
    if not os.path.exists(directory):
        os.mkdir(directory)
# Generate output files for BAM processing
args['bamPrf'] = os.path.join(args['bamDir'], args['name'])
args['alnBam'] = args['bamPrf'] + '.int.bam'
args['alnLog'] = args['bamPrf'] + '.align.log'
args['mdpBam'] = args['bamPrf'] + '.bam'
args['mdpLog1'] = args['bamPrf'] + 'mdup.1.log'
args['mdpLog2'] = args['bamPrf'] + 'mdup.2.log'
# Generate output files for BED processing
args['bedPrf'] = os.path.join(args['bedDir'], args['name'])
args['intBed'] = args['bedPrf'] + '.int.bed'
args['intLog'] = args['bedPrf'] + '.bam2bed.log'
args['srtBed'] = args['bedPrf'] + '.bed'
args['bedGrp'] = args['bedPrf'] + '.bedgraph'
args['bedLog'] = args['bedPrf'] + '.bedprocess.log'
# Generate object to store commands
moabJobs = moab.moabJobs()
# Generate command to perform alignment and add to job dictionary
alignCommand = bowtie2Align(index = args['<bwt2index>'],
    outFile = args['alnBam'], read1 = args['read1'], read2 = args['read2'],
    bowtie2Path = paths['bowtie2'], threads = args['--threads'],
    readGroup = args['name'], sampleName = args['name'],
    libraryID = args['prefix'], discordant = False, mixed = False, upto = 100000,
    maxInsert = 2000, check = True, samtoolsPath = paths['samtools'])
alignCommandID = moabJobs.add(alignCommand, stdout = args['alnLog'],
    stderr = args['alnLog'], processor = args['--threads'])
# Generate command to deduplicate BAM and add to job dictionary
dedupCommand = picard.markDuplicates(inBam = args['alnBam'],
    outBam = args['mdpBam'], logFile = args['mdpLog1'],
    picardPath = paths['picardPath'], javaPath = paths['javaPath'],
    removeDuplicates = False, delete = True)
dedudCommandID = moabJobs.add()
# Generate command to create bed file
bam2bedCommand = '%s %s %s %s --minMapQ %s --size %s--rmDup %s' %(
    paths['python'], paths['bam2bed'], args['mdpBam'], args['intBed'],
    args['--minMapQ'], args['--size'], args['--rmDup'])
# Generate command to sort bed file
sortCommand = bedtools.sortBed(inBed = args['intBed'], outBed = args['srtBed'],
    path = paths['bedtools'], delete = True)
# Generate command to create bedgraph file
bgCommand = bedtools.bed2bedGraph(inBed = args['srtBed'], outBG = args['bedGrp'],
    chrFile = args['<chrFile>'], path = paths['bedtools'], delete = False)
