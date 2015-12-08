from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import random
import multiprocessing
import subprocess
from general_functions import writeFile

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

class readProcessObject(object):
    ''' Creates object to handle reading fastq files in a seperate
    process. Function takes one argument:
    
    1)  fastqFile - Full path to fastq file
    
    Object has two functions:
    
    1)  next - Returns the next fastq read as a string. Raises an
        EOFError if there is no futher reads.
    2)  close - Closes process and the pipe communicating with it.
    
    '''
    
    def __init__(self, fastqFile):
        self.pipes = multiprocessing.Pipe(False)
        self.process = multiprocessing.Process(
            target = readToPipe,
            args = (fastqFile, self.pipes)
        )
        self.process.start()
        self.pipes[1].close()
    
    def next(self):
        return(self.pipes[0].recv())
    
    def close(self):
        self.pipes[0].close()
        self.process.join()

def randomPair(fastqIn1, fastqIn2, fastqOut1, fastqOut2, number):
    ''' Extract random paired end read '''
    # Count reads in file
    if fastqIn1.endswith('.gz'):
        readCount = subprocess.check_output(['zgrep', '-c', '$', fastqIn1])
    else:
        readCount = subprocess.check_output(['grep', '-c', '$', fastqIn1])
    readCount = int(readCount.strip()) / 4
    # Select which reads to extract
    random.seed(1)
    selected = random.sample(
        range(0,readCount),
        min(readCount,number)
    )
    selected.sort(reverse = True)
    # Create processes to read input FASTQ files
    read1In = readProcessObject(fastqIn1)
    read2In = readProcessObject(fastqIn2)
    # Create processes to write output FASTQ files
    read1Out = writeFile.writeProcess(fastqOut1)
    read2Out = writeFile.writeProcess(fastqOut2)
    # Select random reads
    nextRead = selected.pop()
    readCount = 0
    # Extract desired reads
    while True:
        # Extract reads or break loop
        try:
            read1 = read1In.next()
            read2 = read2In.next()
        except EOFError:
            break
        # Find the next random read
        if readCount == nextRead:
            # Check for matching read names and write reads
            read1Name = read1.split('\n',1)[0].split(' ',1)[0]
            read2Name = read2.split('\n',1)[0].split(' ',1)[0]
            if read1Name == read2Name:
                read1Out.add(read1)
                read2Out.add(read2)
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
    read1In.close()
    read2In.close()
    read1Out.close()
    read2Out.close()

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
