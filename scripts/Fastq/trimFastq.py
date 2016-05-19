"""trimFastq.py

Usage:
    
    trimFastq.py <inprefix> <outprefix> [--quality=<quality>]
        [--adapter=<adapter>] [--path=<path>]
    
Options:
    
    --quality=<quality>  Trimming quality [default: 20]
    --adapter=<adapter>  Adapter sequence [default: AGATCGGAAGAGC]
    --path=<path>        Path to cutadapt [default: cutadapt]
    --help               Output this message
    
"""
# Import required modules
import os
from ngs_python.fastq import fastqFind, fastqTrim
from general_python import docopt, toolbox, moab
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
args['--quality'] = int(args['--quality'])
inDir, inPrefix = os.path.split(args['<inprefix>'])
toolbox.checkArg(args['--path'], 'exc')
# Extract fastq files and generate output file names
read1In, read2In = fastqFind.findFastq(prefix = inPrefix, dirList = [inDir],
    pair = True, gzip = True)
read1Out = args['<outprefix>'] + '.R1.fastq.gz'
read2Out = args['<outprefix>'] + '.R2.fastq.gz'
trimLog = args['<outprefix>'] + '.log'
# Generate and submit trim command
trimCommand = fastqTrim.cutadaptTrimPaired(read1In = read1In,
    read2In = read2In, read1Out = read1Out, read2Out = read2Out,
    quality = args['--quality'], adapter = 'AGATCGGAAGAGC', length = 25,
    path = args['--path']
)
jobID = moab.submitJob(trimCommand, stdout = trimLog, stderr = trimLog)
print jobID
