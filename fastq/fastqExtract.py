import gzip
import random
import multiprocessing
import subprocess
import cStringIO
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def readFastqProcess(fastq, pipes, shell = True):
    ''' Function to open FASTQ files, parse entries and add individual
    read to a pipe. Function takes three arguments:
    
    1)  fastq - Full path to FASTQ file.
    2)  pipes - Send enabled pipe.
    3)  shell - Boolean indicating whether to use shell 'zcat' command
        to read FASTQ file.
    
    '''
    # Loop through fastq files and add contents to pipe
    for f in fastqFile:
        # Create fastq handle
        if shell and fastq.endswith('.gz'):
            sp = subprocess.Popen(["zcat", f], stdout = subprocess.PIPE)
            fh = cStringIO.StringIO(sp.communicate()[0])
        elif f.endswith('.gz'):
            fh = gzip.open(f)
        else:
            fh = open(f)
        # Create output and add to pipe
        for read in FastqGeneralIterator(fh):
            pipe.send(read)
        # Add poisin pill to pipe and close
        pipe.close()
        # Try to close the fastq handle
        try:
            fh.close()
        except AttributeError:
            continue

def fastqPairIterator(fastq1, fastq2, shell = True):
    # Create read1 pipes and process
    pipe1Recv, pipe1Send = multiprocessing.Pipe(False)
    process1 = multiprocessing.Process(
        target = readFastqProcess,
        args = (fastq1, pipe1Send, shell) 
    )
    process1.start()
    pipe1Send.close()
    # Create read2 pipes and process
    pipe2Recv, pipe2Send = multiprocessing.Pipe(False)
    process2 = multiprocessing.Process(
        target = readFastqProcess,
        args = (fastq2, pipe2Send, shell) 
    )
    process2.start()
    pipe2Send.close()
    # Create iterator from pipes
    while True:
        try:
            read1 = pipe1Recv.recv()
            read2 = pipe2Recv.recv()
        except EOFError:
            break
        yield((read1, read2))

def fastqCount(fastqFile, shell = True):
    ''' Creates a generator that parses FASTQ files using the 
    FastqGeneralIterator in Bio.SeqIO. Function takes one argument:
    
    1)  fastqFile - Path to FASTQ file or a list of paths.
    
    '''
    # Process input file(s)
    if isinstance(fastqFile, str):
        fastqFile = [fastqFile]
    count = 0
    # Loop through fastq files and add contents to pipe
    for f in fastqFile:
        # Create open file command
        if shell and f.endswith('.gz'):
            sp = subprocess.Popen(["zcat", f], stdout = subprocess.PIPE)
            fh = cStringIO.StringIO(sp.communicate()[0])
        elif f.endswith('.gz'):
            fh = gzip.open(f)
        else:
            fh = open(f)
        # Create output and add to pipe
        for read in FastqGeneralIterator(fh):
            count += 1
        try:
            fh.close()
        except AttributeError:
            continue
    # Return count
    return(count)

def fastqGenerator(fastqFile, shell = True):
    ''' Creates a generator that parses FASTQ files using the 
    FastqGeneralIterator in Bio.SeqIO. Function takes one argument:
    
    1)  fastqFile - Path to FASTQ file or a list of paths.
    
    '''
    # Process input file(s)
    if isinstance(fastqFile, str):
        fastqFile = [fastqFile]
    # Loop through fastq files and add contents to pipe
    for f in fastqFile:
        # Create open file command
        if shell and f.endswith('.gz'):
            sp = subprocess.Popen(["zcat", f], stdout = subprocess.PIPE)
            fh = cStringIO.StringIO(sp.communicate()[0])
        elif f.endswith('.gz'):
            fh = gzip.open(f)
        else:
            fh = open(f)
        # Create output and add to pipe
        for title, seq, qual in FastqGeneralIterator(fh):
            yield('@' + title + '\n' + seq + '\n' '+' '\n' + qual)
        try:
            fh.close()
        except AttributeError:
            continue
    # Raise stop iteration
    raise StopIteration

def readToPipe(fastqFile, pipes, shell = True):
    ''' Function extract reads from a FASTQ file using fastqGenerator
    and adds reads to a pipe. The function closes the pipe when all
    FASTQ reads are processed. Function takes two arguments:
    
    1)  fastqFile - Path to FASTQ file or a list of paths.
    2)  pipes - A pair of multiprocessing connection objects created
        using the multiprocessing.Pipe function.
    
    '''
    # Close unsused receive pipe
    pipes[0].close()
    # Add reads to pipe
    for f in fastqGenerator(fastqFile, shell):
        pipes[1].send(f)
    # Close send pipe
    pipes[1].close()

class readFastq(object):
    ''' Creates object to handle reading fastq files. The FASTQ files
    are parsed using the FastqGeneralIteratorFunction from the Bio.SeqIO
    module. Function takes one argument:
    
    1)  fastqFile - Full path to fastq file
    
    Object has two functions:
    
    1)  next - Returns the next fastq read as a string. Raises an
        EOFError if there is no futher reads.
    2)  close - Closes process and the pipe communicating with it.
    
    '''
    
    
    def __init__(self, fastqFile, shell = True):
        # Store FASTQ files and create generator
        self.fastqFile = fastqFile
        self.shell = shell
        self.generator = fastqGenerator(self.fastqFile, self.shell)
     
    def __iter__(self):
        return(self)
    
    def next(self):
        return(self.generator.next())

class readFastqProcess(object):
    ''' Creates object to handle reading fastq files in a seperate
    process. Function takes one argument:
    
    1)  fastqFile - Full path to fastq file
    
    Object has two functions:
    
    1)  next - Returns the next fastq read as a string. Raises an
        EOFError if there is no futher reads.
    2)  close - Closes process and the pipe communicating with it.
    
    '''
    
    def __init__(self, fastqFile, shell = True):
        # Store suppplied variables
        self.fastqFile = fastqFile
        self.shell = shell
        # Create pipes
        self.pipes = multiprocessing.Pipe(False)
        # Create and start process
        self.process = multiprocessing.Process(
            target = readToPipe,
            args = (fastqFile, self.pipes, self.shell)
        )
        self.process.start()
        # Close input pipe
        self.pipes[1].close()
    
    def __enter__(self):
        return(self)
    
    def __iter__(self):
        return(self)
    
    def next(self):
        try:
            return(self.pipes[0].recv())
        except EOFError:
            raise StopIteration
    
    def close(self):
        self.pipes[0].close()
        self.process.join()
    
    def __del__(self):
        self.close()
    
    def __exit__(self, type, value, traceback):
        self.close()

def randomPair(fastqIn1, fastqIn2, fastqOut1, fastqOut2, number, shell = True):
    ''' Extract random paired end read '''
    # Count reads in file
    readCount = fastqCount(fastqIn1)
    # Select which reads to extract
    random.seed(1)
    selected = random.sample(
        range(0,readCount),
        min(readCount,number)
    )
    selected.sort(reverse = True)
    # Select random reads
    nextRead = selected.pop()
    readCount = 0
    # Create pbjects
    read1In = readFastqProcess(fastqIn1, shell)
    read2In = readFastqProcess(fastqIn2, shell)
    read1Out = writeFile.writeFileProcess(fastqOut1, shell)
    read2Out = writeFile.writeFileProcess(fastqOut2, shell)
    # Loop through fastq files
    for read1, read2 in itertools.izip(read1In, read2In):
        # Find the next random read
        if readCount == nextRead:
            # Check for matching read names and write reads
            read1Name = read1.split('\n',1)[0].split(' ',1)[0]
            read2Name = read2.split('\n',1)[0].split(' ',1)[0]
            if read1Name == read2Name:
                read1Out.add(read1 + '\n')
                read2Out.add(read2 + '\n')
            else:
                raise IOError('Input FASTQ files are not paired')
            # Extract next read number
            try:
                nextRead = selected.pop()
            except IndexError:
                nextRead = 0
        readCount += 1
    # Check that all reads have been extracted
    read1In.close()
    read2In.close()
    read1Out.close()
    read2Out.close()
    if selected:
        raise ValueError('Not all read numbers extracted')

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
