from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import random
import multiprocessing

def readToPipe(fastqFile, pipes):
    ''' Function parses a FASTQ file using Bio.SeqIO and add individual
    fASTQ reads to a supplied multiprocessing pipe. The function closes
    the pipe when all input FASTQ files are processed. Function takes
    two arguments:
    
    1)  fastqFile - A string of the FASTQ file path or a python list
        containing a series of FASTQ file paths.
    2)  pipes - A pair of multiprocessing connection objects created
        using the multiprocessing.Pipe('False') command.
            
    '''
    # Close unsused receive pipe
    pipes[0].close()
    # Process input file(s)
    if isinstance(fastqFile, str):
        fastqFile = [fastqFile]
    # Loop through fastq files and add contents to pipe
    for f in fastqFile:
        # Create open file command
        if f.endswith('.gz'):
            readFile = gzip.open
        else:
            readFile = open
        # Create output and add to pipe
        for title, seq, qual in FastqGeneralIterator(readFile(f)):
            pipes[1].send('@' + title + '\n' + seq + '\n' '+' '\n' + qual)
    # Close send pipe
    pipes[1].close()

def readToPipeProcess(fastqFile):
    ''' Function creates a multiprocessing process and pip. The process
    parses FASTQ files and passes individual entries into the pipe.
    Function returns the process and the end of the pipe from which
    individual FASTQ entries can be extracted. The pipe will return
    EOFError when all entries have been extracted. At this point the
    process can be joined. Function takes one argument:
    
    1)  fastqFile - A string of the FASTQ file path or a python list
    containing a series of FASTQ file paths.
    
    '''
    # Create pipes and process
    pipes = multiprocessing.Pipe(False)
    process = multiprocessing.Process(
        target = readToPipe,
        args = (fastqFile, pipes)
    )
    process.start()
    pipes[1].close()
    # Return pipe
    return(process, pipes[0])

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
    while True:
        try:
            read = pipes[0].recv()
        except EOFError:
            break
        outFile.write(read)
    # Close files and pipes
    outFile.close()
    pipes[0].close()

def writeFromPipeProcess(fastqFile):
    # Create pipes and process
    pipes = multiprocessing.Pipe('False')
    process = multiprocessing.Process(
        target = writeFromPipe,
        args = (fastqFile, pipes)
    )
    process.start()
    pipes[0].close()
    # Return process and pipes
    return(process, pipes[1])

def writeFromPipe2(fileName, pipes):
    # Close unused pipes
    pipes[1].close()
    # Open output file
    fileOut = gzip.open(fileName, 'wb')
    sp = subprocess.Popen('gzip', stdout = fileOut,
        stdin = subprocess.PIPE, bufsize = -1)
    count = 0
    # Write to output file
    while True:
        try:
            read = pipes[0].recv()
        except EOFError:
            break
        sp.stdin.write(read)
        count += 1
        if not (count % 100000):
            print str(count) + ' written'
    # Close files and pipes
    print 'Write loop closed'
    print sp.poll()
    fileOut.close()
    pipes[0].close()
    print "Finished writing"

def randomPair(read1In, read2In, read1Out, read2Out, number):
    ''' Extract random paired end read '''
    # Count reads in file
    if read1In.endswith('.gz'):
        readCount = subprocess.check_output(['zgrep', '-c', '$', read1In])
    else:
        readCount = subprocess.check_output(['grep', '-c', '$', read1In])
    readCount = int(readCount.strip()) / 4
    # Select which reads to extract
    random.seed(1)
    selected = random.sample(
        range(0,readCount),
        min(readCount,number)
    )
    selected.sort(reverse = True)
    # Create processes to read input FASTQ files
    read1InProcess, read1InPipe = readToPipeProcess(read1In)
    read2InProcess, read2InPipe = readToPipeProcess(read2In)
    # Create processes to write output FASTQ files
    read1OutProcess, read1OutPipe = writeFromPipeProcess(read1Out)
    read2OutProcess, read2OutPipe = writeFromPipeProcess(read2Out)
    # Select random reads
    nextRead = selected.pop()
    readCount = 0
    # Extract desired reads
    while True:
        # Extract reads or break loop
        try:
            read1 = read1InPipe.recv()
            read2 = read2InPipe.recv()
        except EOFError:
            break
        # Find the next random read
        if readCount == nextRead:
            # Check for matching read names and write reads
            read1Name = read1.split('\n',1)[0].split(' ',1)[0]
            read2Name = read2.split('\n',1)[0].split(' ',1)[0]
            if read1Name == read2Name:
                read1OutPipe.send(read1)
                read2OutPipe.send(read2)
            else:
                raise IOError('Input FASTQ files are not paired')
            # Extract next read number
            try:
                nextRead = selected.pop()
            except IndexError:
                nextRead = 0
        readCount += 1
    # Check that all reads have been extracted
    if selected:
        raise ValueError('Not all read numbers extracted')
    # Close pipes and processes
    read1InPipe.close()
    read2InPipe.close()
    read1OutPipe.close()
    read2OutPipe.close()
    read1InProcess.join()
    read2InProcess.join()
    read1OutProcess.join()
    read2OutProcess.join()

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
