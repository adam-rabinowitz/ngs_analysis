"""trimFastq.py

Usage:
    
    trimFastq.py <samplename> <inprefix> <outdir> <index>
        [--bowtie2=<bowtie2>] [--rsem=<rsem>] [--threads=<threads>] [--singleend]
        [--forprob=<forprob>]
    
Options:
    
    --bowtie2=<bowtie2>  Path to bowtie2 [default: bowtie2]
    --forprob=<forprob>  Forward probability [default: 0.5]
    --rsem=<rsem>        Path to rsem-calculate-expression [default: rsem-calculate-expression]
    --threads=<threads>  Number of threads [default: 4]
    --help               Output this message
    
"""
# Import required modules
import os
from ngs_python.fastq import fastqFind, fastqAlign
from general_python import docopt, toolbox, moab
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--threads'] = int(args['--threads'])
args['--forprob'] = float(args['--forprob'])
toolbox.checkArg(args['--forprob'], 'num', mn = 0, mx = 1)
inDir, inPrefix = os.path.split(args['<inprefix>'])
outDir = os.path.join(args['<outdir>'], args['<samplename>'])
if not os.path.isdir(outDir):
    os.mkdir(outDir)
outPrefix = os.path.join(outDir, args['<samplename>'])
outLog = outPrefix + '.rsem.log'
# Extract fastq files and generate output file names
read1, read2 = fastqFind.findFastq(prefix = inPrefix, dirList = [inDir],
    pair = True, gzip = True)
rsemCommand = fastqAlign.rsemBowtie2Align(index = args['<index>'],
    outPrefix = outPrefix, read1 = read1, read2 = read2,
    rsemPath = args['--rsem'], bowtie2Path = args['--bowtie2'],
    threads = args['--threads'], forProb = args['--forprob'],
    genomeBam = False, estimateRspd = False, check = True)
jobID = moab.submitJob(rsemCommand, processor = args['--threads'],
    stdout = outLog, stderr = outLog)
print args['<samplename>'], jobID
