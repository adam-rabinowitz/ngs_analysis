import os
import re
import sys
import pysam
import multiprocessing
import subprocess
import collections
from github.general_functions import iohandle

def extract_random(read1In, read2In, read1Out, read2Out, number = 100000):
    ''' This function generates and reutrns a command to extracts random reads
    from paired fastq files. Function takes four arguments:
    
    1)  read1In - Path to read1 gzipped FASTQ input file.
    2)  read2In - Path to read2 gzipped FASTQ input file.
    3)  read1Out - Path to read1 gzipped FASTQ output file.
    4)  read2Out - Path to read2 gzipped FASTQ output file.
    5)  Number of reads to output.
    
    '''
    # Paste files
    step1 = 'paste <(zcat %s) <(zcat %s)' %(
        read1In,
        read2In
    )
    # Merge lines by group of 4
    step2 = 'awk \'{ printf("%s",$0); n++; if(n%4==0) { printf("\\n");} else'+\
        '{ printf("\\t\\t");} }\''
    # Shuffle file and take desired number 
    step3 = 'shuf | head -%s' %(
        number
    )
    # restore delimiter
    step4 = 'sed \'s/\\t\\t/\\n/g\''
    # Print output
    #step5 = 'awk -F \'\\t\' \'{print $1 > "%s"; print $2 > "%s"}\'' %(
    #    read1Out,
    #    read2Out
    #)
    step5 = 'awk -F \'\\t\' \'{print $1 | "gzip > %s"; print $2 | "gzip > %s"}\'' %(
        read1Out,
        read2Out
    )
    # Combine commands
    completeCommand = '%s | %s | %s | %s | %s' %(
        step1,
        step2,
        step3,
        step4,
        step5
    )
    # Return command
    return(completeCommand)


################################################################################
## cutadapt_trim_paired
################################################################################
# Define function
def cutadaptTrimPaired(read1In, read2In, read1Out, read2Out, quality = 20,
    adapter = 'AGATCGGAAGAGC', length = 25, path = 'cutadapt'
):
    ''' This function peforms quality trimming on paired end FASTQ files
    using the cutadapt package. Function is built for version 1.7.1 of
    cutadapt which uses a two stage process to trim FASTQ file. The
    default adapter sequence is common to both ends of the standard
    Illumina adapter. Output file names ending in '.gz' will cause the
    output to be gzipped compressed. Function takes  8 arguments:
    
    1)  Read1 input fastq file.
    2)  Read2 input fastq file.
    3)  Read1 output fastq file.
    4)  Read2 output fastq file.
    5)  Minimum base quality; Default = 20.
    6)  Adapter sequence; Default = 'AGATCGGAAGAGC'.
    7)  Minimum read length; Default = 25. 
    8)  Path to cutadapt; Default = 'cutadapt'.
    
    '''
    # Create temporary FASTQ files
    tempRead1 = read1Out + '.tmp.fastq'
    tempRead2 = read2Out + '.tmp.fastq'
    # Generate trimming commands
    trim1Command = [path, '-q', str(quality), '-a', adapter, '-e',
        '0.1', '--minimum-length', str(length), '-O', '1', '-o', tempRead1,
        '-p', tempRead2, read1In, read2In]
    trim2Command = [path, '-q', str(quality), '-a', adapter, '-e',
        '0.1', '--minimum-length', str(length), '-O', '1', '-o', read2Out,
        '-p', read1Out, tempRead2, tempRead1]
    # Join and return command
    trimCommand = '%s; %s; rm %s %s' %(
        ' '.join(trim1Command),
        ' '.join(trim2Command),
        tempRead1,
        tempRead2
    )
    return(trimCommand)
