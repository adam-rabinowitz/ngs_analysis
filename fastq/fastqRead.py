
import gzip
import random
import multiprocessing
import subprocess
import cStringIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class readFastq(object):

    def __init__(
        self, fastq, shell = True, number = None, sample = None
    ):
        # Store variables
        self.fastq = fastq
        self.shell = shell
        self.number = number
        self.sample = sample
        # Create entries
        if number:
            self.entries = self.genEntries(
                self.number, self.sample
            )
        else:
            self.entries = None
        # Create process
        self.process, self.pipe = self.createProcess(
            self.fastq, self.shell, self.entries
        )
    
    def genEntries(self, number, sample = None):
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
        # Retrun values
        return(entries)
    
    def readProcess(
        self, fastq, pipe, shell = True, entries = None
    ):
        # Process pipe
        pipeRecv, pipeSend = pipe
        pipeRecv.close()
        # Create fastq handle
        if fastq.endswith('.gz'):
            fh = gzip.open(fastq)
        else:
            fh = open(fastq)
        # If number specified extract desired reads
        if entries:
            # Create count and sort numbers
            count = 0
            nextRead = entries.pop()
            # Loop through all reads and send match to pipe
            for read in FastqGeneralIterator(fh):
                if count == nextRead:
                    pipeSend.send(read)
                    try:
                        nextRead = entries.pop()
                    except IndexError:
                        break
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
            self.entries.reverse()
    
    def createProcess(
        self, fastq, shell, entries
    ):
        ''' Funtion to generate process to read FASTQ file '''
        # Create pipe and process
        pipe = multiprocessing.Pipe(False)
        process = multiprocessing.Process(
            target = self.readProcess,
            args = (fastq, pipe, shell, entries)
        )
        process.start()
        pipe[1].close()
        # Return pipe and process
        return(process, pipe)
    
    def __iter__(self):
        return(self)
    
    def next(self):
        ''' Returns next element in pipe or raises StopIteration error
        when pipe is empty
        '''
        try:
            return(self.pipe[0].recv())
        except EOFError:
            raise StopIteration

x = readFastq(
    '/farm/scratch/rs-bio-lif/rabino01/Elza/fastqFiles/NGS-10251_0611_L001_R1.fastq.gz',
    number = 10,
    sample = 100
)
print x
for v in x:
    print v
