from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

def readToPipe(fastqFile, pipes):
    ''' Function parses a FASTQ file using Bio.SeqIO and creates a list
    object where each element is a line of a single FASTQ entry. The
    list elements are sent down the pipe. Function takes two arguments:
    
    1)  fastqFile - Full path to input FASTQ file.
    2)  pipes - Tuple containg paired multiprocessing.Connection
        objects. Pair should be unidirectional with first used to
        receive messages snd the second to send messages.
            
    '''
    # Close unsused pipes
    pipes[0].close()
    # Create command to open file
    if fastqFile.endswith('.gz'):
        openFile = gzip.open
    else:
        openFile = open
    # Create output
    for title, seq, qual in FastqGeneralIterator(openFile(fastqFile)):
        pipes[1].send(['@' + title, seq, '+', qual])
    pipes[1].close()

def writeFromPipe(fileName, pipes):
    ''' Write FASTQ or FASTA file 
    1)  fastqFile - Full path to output FASTQ file
    '''
    # Close unused pipes
    pipes[1].close()
    # Open outfile dependant on prefix
    if fileName.endswith('.gz'):
        outFile = gzip.open(fileName, 'w')
    else:
        outFile = open(fileName, 'w')
    # Write to output file
    count = 0
    while True:
        try:
            read = pipes[0].recv()
        except EOFError:
            break
        outFile.write('\n'.join(read) + '\n')
        count += 1
        if not (count % 100000):
            print str(count) + ' written'
    # Close files and pipes
    print 'Write loop closed'
    outFile.close()
    pipes[0].close()
    print "Finished writing"

def writeFromPipe2(fileName, pipes):
    # Close unused pipes
    pipes[1].close()
    # Open output file
    fileOut = gzip.open(fileName, 'wb')
    sp = subprocess.Popen('gzip', stdout = fileOut,
        stdin = subprocess.PIPE)
    count = 0
    # Write to output file
    while True:
        try:
            read = pipes[0].recv()
        except EOFError:
            break
        sp.stdin.write('\n'.join(read) + '\n')
        count += 1
        if not (count % 100000):
            print str(count) + ' written'
    # Close files and pipes
    print 'Write loop closed'
    sp.terminate()
    fileOut.close()
    pipes[0].close()
    print "Finished writing"

import os
import subprocess
from cStringIO import StringIO

def read2Pipe2(fastqFile, pipe):
    p = subprocess.Popen(["zcat", fastqFile], stdout = subprocess.PIPE)
    f = StringIO(p.communicate()[0])
    assert p.returncode == 0
    for title, seq, qual in FastqGeneralIterator(f):
        pipe.send(['@' + title, seq, '+', qual])
    pipe.close()
    

def extractRandom(read1In, read2In, read1Out, read2Out, number = 100000):
    ''' This function generates and reutrns a command to extracts random
    reads from paired fastq files. Function takes five arguments:
    
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
    # Shuffle file
    step3 = 'shuf'
    # Extract desired number of reads
    step4 = 'head -%s' %(
        number
    )
    # restore delimiter
    step5 = 'sed \'s/\\t\\t/\\n/g\''
    # Print output
    step6 = 'awk -F \'\\t\' \'{print $1 | "%s"; print $2 | "%s"}\'' %(
        'gzip > ' + read1Out,
        'gzip > ' + read2Out
    )
    # Combine and return command
    completeCommand = '%s | %s | %s | %s | %s | %s' %(
        step1, step2, step3, step4, step5, step6)
    return(completeCommand)
