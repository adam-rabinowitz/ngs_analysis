import gzip
import multiprocessing
import subprocess
import os
import random
from Bio.SeqIO.QualityIO import FastqGeneralIterator

class readFastq(object):
    ''' Class functions as an iterator to extract reads from single or paired
    FASTQ files. Classes uses multiprocessing to speed up the extraction of
    the FASTQ entries and the FastqGeneralIterator from Bio.SeqIO to parse
    indivdual reads.
    '''
    
    def __init__(
        self, fastq1, fastq2 = None, shell = True, number = None, sample = None
    ):
        ''' Function to initialise object. Function takes five arguments:
        
        1)  fastq1 - Full path to FASTQ file.
        2)  fastq2 - Full path to paired FASTQ file (optional).
        3)  shell - Boolean indicating whether to use shell zcat command
            to read gzipped input.
        4)  number - Number of reads to extract.
        5)  sample - Number of reads from which to sample the desired number
            of read. Only used when 'number' also set.
        
        '''
        # Store fastq files
        self.fastq_list = [fastq1]
        if fastq2:
            self.pair = True
            self.fastq_list.append(fastq2)
        else:
            self.pair = False
        # Store additional variables
        self.shell = shell
        self.number = number
        self.sample = sample
        # Create entries
        if number:
            self.entries = self.__gen_entries()
        else:
            self.entries = None
        # Create list to store process data
        self.process_list = []
    
    def __gen_entries(self):
        ''' Function utilises the self.number and self.sample values to
        generate the index of the desired entries within the FASTQ file.
        '''
        # Generate random entries from provided sample
        if self.sample:
            # Check values
            if self.number > self.sample:
                raise ValueError('Number must be smaller than sample')
            # Create random sample
            random.seed(1)
            entries = random.sample(range(self.sample), self.number)
            entries.sort(reverse = True)
        # Or select first desired entries
        else:
            entries = range(number)
            entries.reverse()
        # Retrun values
        return(entries)
    
    def __read_process(self, fastq, pend):
        ''' Function to generate a process to read FASTQ files. Extracted reads
        will be sent doen the supplied multiprocessing pipe. If self.shell is
        True then gzipped input files will be read using the zcat command in the
        shell. If self.entries is set then only the FASTQ entries at the desired
        indices will be sent down the pipe. FASTQ entries are generated as three
        element tuples consisting of: read name, base calls, base qualities.
        Function takes two arguments:
        
        1)  fastq - Full path to the FASTQ file to read
        2)  pend - End of multiprocessing pipe down which reads will be sent.
        '''
        # Create fastq handle
        if self.shell and fastq.endswith('.gz'):
            sp = subprocess.Popen(['zcat', fastq], stdout = subprocess.PIPE,
                bufsize = 1)
            fh = sp.stdout
        elif fastq.endswith('.gz'):
            fh = gzip.open(fastq)
        else:
            fh = open(fastq)
        # If number specified extract desired reads
        if self.entries:
            # Create count and sort numbers
            count = 0
            nextRead = self.entries.pop()
            # Loop through all reads and send match to pipe
            for read in FastqGeneralIterator(fh):
                if count == nextRead:
                    # Send data or break iteration
                    try:
                        self.__send(read, pend)
                    except StopIteration:
                        entries = []
                    # Extract next read or break iteration
                    try:
                        nextRead = self.entries.pop()
                    except IndexError:
                        break
                count += 1
        # Else extract all reads
        else:
            for read in FastqGeneralIterator(fh):
                # Send data or break iteration
                try:
                    self.__send(read, pend)
                except StopIteration:
                    break
        # Close subprocess and file handle
        if self.shell and fastq.endswith('.gz'):
            sp.terminate()
        fh.close()
        # Raise error if not all desired reads extracted
        if self.entries:
            raise IOError('Not all reads extracted from %s' %(fastq))
            self.entries.reverse()
    
    def __send(self, read, pend):
        # Check pipe for incoming signal
        if pend.poll():
            # Extract signal and process
            recv = pend.recv()
            if recv is None:
                raise StopIteration('Termination signal received')
            else:
                raise ValueError('Unknown signal received')
        # Else add read to pipe
        else:
            pend.send(read)
    
    def start(self):
        ''' Function creates processes to read the FASTQ files listed in
        self.fastq_list using self.__read_process. Function creates two
        element tuples consisting of the process and pipe end with which
        to recevie FASTQ reads from the process. The tuples are stored in
        self.process_list.
        '''
        # Close active processes
        self.close()
        # Loop through fastq files:
        for fastq in self.fastq_list:
            # Create pipe and process
            pend_child, pend_parent = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__read_process,
                args = (fastq, pend_child)
            )
            process.start()
            pend_child.close()
            # Store process data
            self.process_list.append((process, pend_parent))
    
    def close(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. The list 
        self.process_list is then emptied.
        '''
        # Loop through processes and terminate
        for process, pend in self.process_list:
            # Add termination signal to pipes
            try:
                pend.send(None)
            except IOError:
                pass
            # Join process and close pipes
            process.join()
            pend.close()
        # Empty process list
        self.process_list = []
    
    def __iter__(self):
        ''' Returns self for iteration '''
        return(self)
    
    def next(self):
        ''' Returns next element in pipe or raises StopIteration '''
        # Return fastq pairs
        if self.pair:
            try:
                return((
                    self.process_list[0][1].recv(),
                    self.process_list[1][1].recv()
                ))
            except EOFError:
                raise StopIteration
        # Else return single reads
        else:
            try:
                return(self.process_list[0][1].recv())
            except EOFError:
                raise StopIteration
    
    def __enter__(self):
        ''' Starts processes at start of with scope '''
        # Close active processes
        self.close()
        # Start new processes
        self.start()
        return(self)
    
    def __exit__(self, exception_type, exception_value, traceback):
        ''' Closes processes at end of with scope '''
        self.close()


class writeFastq(object):
    ''' An object that uses multiprocessing processes to parralelize the
    writing of FASTQ files.
    '''
    
    def __init__(self, fastq1, fastq2 = None, shell = True):
        ''' Function to initialise object. Function takes three arguments:
        
        1)  fastq1 - Full path to FASTQ file.
        2)  fastq2 - Full path to paired FASTQ file (optional).
        3)  shell - Boolean indicating whether to use shell gzip command
            to write gzipped output.
        
        '''
        # Store fastq files
        self.fastq_list = [fastq1]
        if fastq2:
            self.fastq_list.append(fastq2)
            self.pair = True
        else:
            self.pair = False
        # Store shell argument
        self.shell = shell
        # Create lists to store pipes and processes
        self.process_list = []
    
    def __write_process(self, fastq, pend):
        ''' Function to generate a process to write FASTQ files. FASTQ reads
        received from the multiprocessing pipe will be written to file. Receipt
        of None will cause the termination of the process. If self.shell is True
        then gzipped output files will be written using the gzip command in the
        shell. Funcion takes two arguments:
        
        1)  fastq - Full path to the FASTQ file to be created
        2)  pend - End of multiprocessing pipe from which reads will be
            extracted.
        '''
        # Write gzip file using shell
        if self.shell and fastq.endswith('.gz'):
            # Create process
            command = 'gzip -c > %s' %(fastq)
            sp = subprocess.Popen(command, shell=True,
                stdin = subprocess.PIPE)
            fh = sp.stdin
        # Write file using python
        elif fastq.endswith('.gz'):
            fh = gzip.open(fastq, 'w')
        else:
            fh = open(fastq, 'w')
        # Extract reads from pipe and write to file
        while True:
            # Extract read
            read = pend.recv()
            # Break loop if read is none
            if read is None:
                break
            # Write read to file
            readString = '{}\n{}\n+\n{}\n'.format(read[0], read[1], read[2])
            fh.write(readString)
        # Close files, pipes and subprocesses
        fh.close()
        if self.shell and fastq.endswith('.gz'):
            sp.communicate()
    
    def start(self):
        ''' Function creates processes to write the FASTQ files listed in
        self.fastq_list using self.__write_process. Function creates two
        element tuples consisting of the process and pipe end with which
        to communicate with the process. The tuples are stored in
        self.process_list.
        '''
        # Close active processes
        self.close()
        # Loop through fastq files
        for fastq in self.fastq_list:
            # Create pipe and process
            pend_child, pend_parent = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = self.__write_process,
                args = (fastq, pend_child)
            )
            process.start()
            pend_child.close()
            # Store pipe end and processes
            self.process_list.append((process, pend_parent))
    
    def close(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. Thelist 
        self.process_list is then emptied.
        '''
        # Extract process and pipes
        for process, pend in self.process_list:
            # Add poisin pill, join process and close pipe
            pend.send(None)
            process.join()
            pend.close()
        # Clear pipe and process list
        self.process_list = []
    
    def __pipe_send(self, read, pend):
        ''' Function to send a FASTQ read down a multiprocessing pipe. If
        the read does not correspond to the expected format alla active
        processes are terminated and an IOError is raised. Function takes
        two arguments:
        
        1)  read - FASTQ read. This should consist of a tuple/list of
            three elements: read name, sequence, quality string.
        2)  pend - End of multiprocessing pipe down which read should be
            sent.
        '''
        # Add read1 to pipe
        if isinstance(read, (tuple,list)) and len(read) == 3:
            pend.send(read)
        else:
            self.close()
            raise IOError('Read must be a list/tuple of three elements')
    
    def write(self, reads):
        ''' Function to write paired reads to paired FASTQ files. Function
        takes one argument:
        
        1)  reads - For single FASTQ files a single read should be provided
            while for paired FASTQ files a tuple/list of 2 reads should be
            provided. Each read should consist of a tuple/list of three
            strings: read name, sequence and quality
        '''
        # Process pairs
        if self.pair:
            # Check input consists of two elements
            if len(reads) != 2:
                self.close()
                raise IOError('Write requires two elements for paired FASTQ')
            # Send elements to write process
            self.__pipe_send(reads[0], self.process_list[0][1])
            self.__pipe_send(reads[1], self.process_list[1][1])
        # Process single reads
        else:
            # Send read to pipe
            self.__pipe_send(reads, self.process_list[0][1])
    
    def __enter__(self):
        ''' Start processes upon entry into with scope '''
        self.start()
        return(self)
    
    def __exit__(self, type, value, traceback):
        ''' Terminate processes upon exit of with scope '''
        self.close()


def random_fastq(number, sample, fastq1In, fastq1Out, fastq2In = None,
    fastq2Out = None, shell = True
):
    ''' Function to extract random reads from a fastq file '''
    with readFastq(fastq1 = fastq1In, fastq2 = fastq2In, number = number,
        sample = sample, shell = shell) as fastqIn:
        with writeFastq(fastq1 = fastq1Out, fastq2 = fastq2Out,
            shell = shell) as fastqOut:
            for read in fastqIn:
                fastqOut.write(read)

random_fastq(
    number = 1000000,
    sample = 60000000,
    fastq1In = '/farm/scratch/rs-bio-lif/rabino01/Elza/fastqFiles/NGS-10251_0611_L001_R1.fastq.gz',
    fastq2In = '/farm/scratch/rs-bio-lif/rabino01/Elza/fastqFiles/NGS-10251_0611_L001_R2.fastq.gz',
    fastq1Out = '/farm/home/rabino01/testOut.R1.fastq.gz',
    fastq2Out = '/farm/home/rabino01/testOut.R2.fastq.gz'
)

