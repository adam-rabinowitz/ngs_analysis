import gzip
import random
import multiprocessing
import subprocess
import cStringIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def readFastqProcess(fastq, pipes, shell = True, entries = None):
    ''' Function to open FASTQ files, parse entries and add individual
    read to a pipe. Function takes three arguments:
    
    1)  fastq - Full path to FASTQ file.
    2)  pipes - Send enabled pipe.
    3)  shell - Boolean indicating whether to use shell 'zcat' command
        to read FASTQ file.
    
    '''
    # Unpack and process pipes
    print 'read process', fastq, entries
    pipeRecv, pipeSend = pipes
    pipeRecv.close()
    # Create fastq handle
    if shell and fastq.endswith('.gz'):
        print 'opening with shell'
        sp = subprocess.Popen(['zcat', fastq], stdout = subprocess.PIPE)
        print 'create subprocess'
        fh = cStringIO.StringIO(sp.communicate()[0])
        print 'create filehandle'
    elif fastq.endswith('.gz'):
        fh = gzip.open(f)
    else:
        fh = open(f)
    print 'fh open'
    # If number specified extract desired reads
    if entries:
        # Create count and sort numbers
        count = 0
        nextRead = entries.pop()
        # Loop through all reads
        for read in FastqGeneralIterator(fh):
            # Find matching count
            if count == nextRead:
                # Send read and extract next number to count
                pipeSend.send(read)
                try:
                    nextRead = entries.pop()
                except IndexError:
                    break
            # Increment read counter
            count += 1
    # Else extract all reads
    else:
        for read in FastqGeneralIterator(fh):
            pipeSend.send(read)
    # Close pipe and file handle
    pipeSend.close()
    fh.close()
    # Raise error if not all desired reads extracted
    if entries:
        raise IOError('Not all reads extracted from %s' %(fastq))

def fastqIterator(fastq1, fastq2 = None, shell = True, number = None,
        sample = None
    ):
    # Generate index of entries required
    if number:
        # Generate random entries from provided sample
        if sample:
            # Check values
            if number > sample:
                raise ValueError('Number must be smaller than sample')
            # Create random sample
            random.seed(1)
            entries = random.sample(range(sample), number)
            entries.sort(reverse = True)
        # Or select first desired entries
        else:
            entries = range(number)
            entries.reverse()
    # Or return no entries
    else:
        entries = None
    # Create read1 pipes and process
    pipes1 = multiprocessing.Pipe(False)
    process1 = multiprocessing.Process(
        target = readFastqProcess,
        args = (fastq1, pipes1, shell, entries) 
    )
    process1.start()
    pipes1[1].close()
    # Create read2 pipes and process
    if fastq2:
        pipes2 = multiprocessing.Pipe(False)
        process2 = multiprocessing.Process(
            target = readFastqProcess,
            args = (fastq2, pipes2, shell, entries) 
        )
        process2.start()
        pipes2[1].close()
    # Create iterator for pairs
    if fastq2:
        while True:
            # Extract values
            try:
                read1 = pipes1[0].recv()
                read2 = pipes2[0].recv()
            # Clean up at end
            except EOFError:
                #pipes1[0].close()
                #pipes2[0].close()
                #process1.join()
                #process2.join()
                break
            # Yield values
            print read1, read2
            yield((read1, read2))
    # Create iterator for single reads
    else:
        while True:
            # Extract values
            try:
                read1 = pipes1[0].recv()
            # Clean up at end
            except EOFError:
                #pipes1[0].close()
                #process1.join()
                break
            # Yield values
            yield(read1)

def writeFastqProcess(fastq, pipes, shell = True):
    ''' Function to write FASTQ files using read passed from a pipe. It is
    expected that the pipe will deliver a list/tuple of three elements of
    the first, third and fourth line of each FASTQ entry. A '+' symbol will
    be inserted into the third line. Function takes three arguments:
    
    1)  fastq - Full path to FASTQ file.
    2)  pipes - Both ends of a pipe from which FASTQ data will be extracted.
    3)  shell - Boolean indicating whether to use shell 'zcat' command
        to read FASTQ file.
    
    '''
    # Unpack pipes and close send end
    pipeRecv, pipeSend = pipes
    pipeSend.close()
    # Write gzip file using shell
    if shell and fastq.endswith('.gz'):
        # Create command
        command = 'gzip -c > %s' %(fastq)
        # Create process
        fh = subprocess.Popen(command, shell=True,
            stdin = subprocess.PIPE)
        # Extract reads from pipe and write to file
        while True:
            # Extract read
            read = pipeRecv.recv()
            # Break loop if none value received
            if read is None:
                break
            # Loop through read:
            readString = '%s\n%s\n+\n%s\n' %(read[0], read[1], read[2])
            fh.stdin.write(readString)
        # Close process
        fh.communicate()
    # Write outfile using python
    else:
        # Open output file
        if fastq.endswith('.gz'):
            fh = gzip.open(fastq, 'w')
        else:
            fh = open(fastq, 'w')
        # Extract reads from pipe and write to file
        while True:
            # Extract read
            read = pipeRecv.recv()
            # Break loop if read is none
            if read is None:
                pipeRecv.close()
                break
            # Loop through read:
            readString = '%s\n%s\n+\n%s\n' %(read[0], read[1], read[2])
            fh.write(readString)
        # Close the file or process
        fh.close()

class writeFastq(object):
    
    def __init__(self, fastq1, fastq2 = None, shell = True):
        # Store fastq files
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        # Create read1 pipes and process
        self.pipes1 = multiprocessing.Pipe(False)
        self.process1 = multiprocessing.Process(
            target = writeFastqProcess,
            args = (self.fastq1, self.pipes1, shell)
        )
        self.process1.start()
        self.pipes1[0].close()
        # Create read2 pipes and process
        if self.fastq2:
            self.pipes2 = multiprocessing.Pipe(False)
            self.process2 = multiprocessing.Process(
                target = writeFastqProcess,
                args = (self.fastq2, self.pipes2, shell)
            )
            self.process2.start()
            self.pipes2[0].close()
        # Set write function
        if self.fastq2:
            self.write = self.writePair
        else:
            self.write = self.writeSingle
    
    def __enter__(self):
        return(self)

    def writeSingle(self, read1):
        # Add reads to pipe
        self.pipes1[1].send(read1)
        
    def writePair(self, read1, read2):
        # Add reads to pipe
        self.pipes1[1].send(read1)
        self.pipes2[1].send(read2)
    
    def close(self):
        # Close pipes and processes
        self.pipes1[1].send(None)
        self.pipes1[1].close()
        self.process1.join()
        if self.fastq2:
            self.pipes2[1].send(None)
            self.pipes2[1].close()
            self.process2.join()
    
    def __exit__(self, type, value, traceback):
        self.close()

fastqOut = writeFastq(
    '/farm/home/rabino01/testOut.R1.fastq.gz',
    '/farm/home/rabino01/testOut.R1.fastq.gz',
    shell = True
)
fastqOut.close()
for x in fastqIterator(
        fastq1 = '/farm/home/rabino01/testIn.R1.fastq.gz',
        fastq2 = '/farm/home/rabino01/testIn.R2.fastq.gz',
        shell = True,
        number = 10
    ):
    print x
    fastqOut.write(x[0], x[1])
fastqOut.close()

## Run script
#if __name__ == "__main__":
#    # Extract and parse argument
#    import argparse
#    parser = argparse.ArgumentParser()
#    parser.add_argument('-I1', type = str, required = True,
#        help = 'Input read1 FASTQ')
#    parser.add_argument('-I2', type = str, help = 'Input read1 FASTQ')
#    parser.add_argument('-O1', type = str, required = True,
#        help = 'Output read1 FASTQ')
#    parser.add_argument('-O2', type = str, help = 'Output read1 FASTQ')
#    parser.add_argument('-N', type = int, required = True,
#        help = 'No. of reads to extract')
#    parser.add_argument('-S', type = int, required = True,
#        help = 'No. of reads to sample')
#    args = parser.parse_args()
#    print args
#    # Extract pairs
#    if args.I2 and args.O2:
#        print 'Opening'
#        outFastq = writeFastq(
#            fastq1 = args.O1,
#            fastq2 = args.O2,
#            shell = True
#        )
#        print 'Writing'
#        for read1, read2 in fastqIterator(
#            fastq1 = args.I1,
#            fastq2 = args.I2,
#            shell = True,
#            number = args.N,
#            sample = args.S
#        ):
#            print read1
#            outFastq.write(read1, read2)
#        outFastq.close()
#    # Raise IOError for invalid entries
#    elif args.I2:
#        raise IOError('No read2 output file provided')
#    elif args.O2:
#        raise IOError('No read2 input file provided')
#    # Extract singletons
#    else:
#        # Sample and write reads
#        with writeFastq(
#            fastq1 = args.O1,
#            shell = True
#        ) as outFastq:
#            for read1, read2 in fastqIterator(
#                fastq1 = args.I1,
#                shell = True,
#                number = args.N,
#                sample = args.S
#            ):
#                print read1
#                outFastq.write(read1, read2)
