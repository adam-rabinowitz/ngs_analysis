import gzip
import multiprocessing
import subprocess
import os
import random

def FastqGeneralIterator(handle): 
    ''' Function taken from Bio.SeqIO.Quality.IO '''
    # We need to call handle.readline() at least four times per record, 
    # so we'll save a property look up each time: 
    handle_readline = handle.readline 
    
    # Skip any text before the first record (e.g. blank lines, comments?) 
    while True: 
        line = handle_readline() 
        if not line: 
            return
        if line[0] == "@": 
           break 
        if isinstance(line[0], int): 
             raise ValueError("Is this handle in binary mode not text mode?") 
   
    while line: 
        if line[0] != "@": 
            raise ValueError( 
                "Records in Fastq files should start with '@' character") 
        title_line = line[1:].rstrip()
        # Will now be at least one line of quality data - in most FASTQ files 
        # just one line! We therefore use string concatenation (if needed) 
        # rather using than the "".join(...) trick just in case it is multiline: 
        seq_string = handle_readline().rstrip() 
        # There may now be more sequence lines, or the "+" quality marker line: 
        while True: 
            line = handle_readline() 
            if not line: 
                raise ValueError("End of file without quality information.") 
            if line[0] == "+": 
                # The title here is optional, but if present must match! 
                second_title = line[1:].rstrip() 
                if second_title and second_title != title_line: 
                    raise ValueError("Sequence and quality captions differ.") 
                break 
            seq_string += line.rstrip()  # removes trailing newlines 
        seq_len = len(seq_string) 

        # Will now be at least one line of quality data... 
        quality_string = handle_readline().rstrip() 
        # There may now be more quality data, or another sequence, or EOF 
        while True: 
            line = handle_readline() 
            if not line: 
                break
            if line[0] == "@": 
                # This COULD be the start of a new sequence. However, it MAY just 
                # be a line of quality data which starts with a "@" character.  We 
                # should be able to check this by looking at the sequence length 
                # and the amount of quality data found so far. 
                if len(quality_string) >= seq_len: 
                    # We expect it to be equal if this is the start of a new record. 
                    # If the quality data is longer, we'll raise an error below. 
                    break 
                    # Continue - its just some (more) quality data. 
                quality_string += line.rstrip() 
        if seq_len != len(quality_string): 
            raise ValueError("Lengths of sequence and quality values differs " 
                " for %s (%i and %i)." 
                % (title_line, seq_len, len(quality_string))) 

        # Return the record and then continue... 
        yield (title_line, seq_string, quality_string) 
    raise StopIteration 

class parseFastq(object):
    
    ''' Class functions as an iterator to extract reads from single or paired
    FASTQ files. Classes uses multiprocessing to speed up the extraction of
    the FASTQ entries and the FastqGeneralIterator from Bio.SeqIO to parse
    indivdual reads.
    '''
    
    def __init__(
        self, fastq1, fastq2 = None, shell = True
    ):
        ''' Function to initialise readFastq object. Checks FASTQ files
        exist and creates lists to store read and write processes.
        
        Args:
            fastq1 (str)- Full path to fastq file.
            fastq2 (str)- Full path to paired fastq file.
            shell (bool)- Whether to use shell to read gzip files.
        
        '''
        # Store fastq files
        self.fastq_list = [fastq1]
        if not fastq2 is None:
            self.pair = True
            self.fastq_list.append(fastq2)
        else:
            self.pair = False
        # Check if fastq files exist
        for fastq in self.fastq_list:
            if not os.path.isfile(fastq):
                raise IOError('File {} could not be found'.format(fastq))
        # Store and create additional variables
        self.shell = shell
        self.read_processes = []
    
    def __read_handle_create(self, fastq):
        ''' Function to create filehandle and subprocess for reading of
        FASTQ files.
        
        Args:
            fastq (str)- Full path to fastq file.
        
        Returns:
            fh - File handle for fastq file.
            sp - zcat subprocess reading the FASTQ file if self.shell=True
                and files ends with '.gz', else None.
        
        '''
        # Create fastq handle and subprocess
        if self.shell and fastq.endswith('.gz'):
            sp = subprocess.Popen(['zcat', fastq], stdout = subprocess.PIPE,
                bufsize = 1)
            fh = sp.stdout
        elif fastq.endswith('.gz'):
            sp = None
            fh = gzip.open(fastq)
        else:
            sp = None
            fh = open(fastq)
        # Return handle and subprocess
        return(fh, sp)
    
    def __send(self, data, conn):
        ''' Function to send data down the end of multiprocessing pipe.
        First checks that termination signal, None, has not been received
        first.
        
        Args:
            data - Data to send down pipe.
            conn - End of pipe.
            
        Raises:
            StopIteration: If termination signal of None is received.
            ValueError: If anythind other than None is received.
       
        '''
        # Check pipe for incoming signal
        if conn.poll():
            # Extract signal and process
            recv = pend.recv()
            if recv is None:
                raise StopIteration('Termination signal received')
            else:
                raise ValueError('Unknown signal received')
        # Add read to connection
        else:
            conn.send(data)
    
    def __read_processes_recv(self):
        ''' Function to receive data from running processes.
        
        Returns:
            data - A list of data received from each process.
        
        Raises:
            IOError - If not all processes return data.
            EOFError - If all processes have no data
        
        '''
        # Extract signal and process
        data = []
        if self.pair:
            try:
                data.append(self.read_processes[0][1].recv())
            except EOFError:
                try:
                    self.read_processes[1][1].recv()
                    raise IOError('Differing number of data points')
                except EOFError:
                    raise EOFError
            try:
                data.append(self.read_processes[1][1].recv())
            except EOFError:
                raise IOError('Differing number of data points')
        else:
            data = data.append(self.read_processes[0][1].recv())
        return(data)
    
    def __read_process_stop(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. The list 
        self.process_list is then emptied.
        '''
        # Loop through processes and terminate
        for process, conn in self.read_processes:
            # Add termination signal to pipes
            try:
                conn.send(None)
            except IOError:
                pass
            # Join process and close pipes
            process.join()
            conn.close()
        # Empty process list
        self.read_processes = []
    
    def __sample_read_process(self, fastq, entries, conn):
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
        # Process pipes
        conn1, conn2 = conn
        conn1.close()
        # Create file handle and subprocess
        fh, sp = self.__read_handle_create(fastq)
        # Create count and sort numbers
        count = 0
        nextRead = entries.pop()
        # Loop through all reads and send match to pipe
        for read in FastqGeneralIterator(fh):
            if count == nextRead:
                # Send data or break iteration
                try:
                    self.__send(read, conn2)
                except StopIteration:
                    break
                # Extract next read or break iteration
                try:
                    nextRead = entries.pop()
                except IndexError:
                    break
            count += 1
        # Clean up
        if sp:
            sp.terminate()
        fh.close()
        conn2.close()
    
    def sample_reads(
            self, number, outFastq1, outFastq2=None, sample=None, seed=1234
        ):
        ''' Function to sample reads from fastq files and write to output
        fastq file.
        
        Args:
            outFastq1 (str)- Full path to output fastq file.
            outFastq2 (str)- Full path to paired output fastq file.
            number (int)- Number of reads to extract.
            sample (int)- Number of reads from which to sample reads.
            seed (int)- Seed for random number generator
        
        '''
        # Check output fastq files
        if self.pair:
            if outFastq2 is None:
                raise ValueError('No output file for 2nd fastq input')
        else:
            if not outFastq2 is None:
                raise ValueError('No 2nd output file required')
        # Create random sample or define index of desired entries
        if sample:
            if number > sample:
                raise ValueError('Number must be smaller than sample')
            random.seed(seed)
            entries = random.sample(range(sample), number)
            entries.sort(reverse = True)
        else:
            entries = range(number)
            entries.reverse()
        # Create pipes and processes to extract fastq and store
        for fastq in self.fastq_list:
            conn = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__sample_read_process,
                args = (fastq, entries, conn)
            )
            process.start()
            conn[1].close()
            self.read_processes.append((process, conn[0]))
        # Create processes to write output reads
        count = 0
        with writeFastq(outFastq1, outFastq2, self.shell) as fastqOut:
            while True:
                try:
                    data = self.__read_processes_recv()
                except EOFError:
                    break
                count += 1
                fastqOut.write(data)
        # Clean up
        self.__read_process_stop()
        # Raise error if desired number of reads not extracted
        if count != number:
            raise IOError('Not all desired reads extracted')
    
    def __check_name_process(self, fastq, conn):
        ''' Function to generate a process to read FASTQ files and extract read
        name and number.
        
        Args:
            fastq (str)- Full path to the FASTQ file to read
            pend - End of multiprocessing pipe down which data will be sent.
        
        Sends:
            name (str)- Read name
            number (int)- Read number
        
        '''
        # Process pipes
        conn1, conn2 = conn
        conn1.close()
        # Loop through fastq files
        fh, sp = self.__read_handle_create(fastq)
        for read in FastqGeneralIterator(fh):
            # Extract read
            name, other = read[0].split(None, 1)
            number = other.split(':', 1)[0]
            # Send data or break iteration
            try:
                self.__send((name, number), conn2)
            except StopIteration:
                break
        # Clean up
        if sp:
            sp.terminate()
        fh.close()
        conn2.close()
    
    def check_names(self):
        ''' Arguments check fastq files and returns count of reads.
        For paired fastq files function checks read names match and
        the read number for fastq files in self.fastq_list are 1 and
        2, respectively. For single fastq file function checks that the
        read number is 1.
        
        Returns:
            count (int) - Number of reads in each FASTQ file.
        
        '''
        # Loop through fastq files:
        for fastq in self.fastq_list:
            # Create pipe and process
            conn = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__check_name_process,
                args = (fastq, conn)
            )
            process.start()
            conn[1].close()
        # Extract data and count reads
        count = 0
        while True:
            try:
                data = conn[0].recv()
            except EOFError:
                break
            count += 1
            # Check names and read number for paired fastq
            if self.pair:
                if (data[0][0] != data [1][0]
                    or data[0][1] != '1'
                    or data[1][1] != '2'):
                    self.__read_process_stop()
                    raise ValueError('Read {} name error'.format(count))
            # Check read number for single fastq
            else:
                if data[1] != '1':
                    self.__read_process_stop()
                    raise ValueError('Read {} name error'.format(count))
        # Stop read processes and return count
        self.__read_process_stop()
        return(count)
    
    def __interleave_read_process(self, fastq, conn):
        ''' Function to generate a process to read FASTQ files and extract read
        name and number.
        
        Args:
            fastq (str)- Full path to the FASTQ file to read
            pend - End of multiprocessing pipe down which data will be sent.
        
        Sends:
            name (str)- Read name
            number (int)- Read number
        
        '''
        # Process pipes
        conn1, conn2 = conn
        conn1.close()
        # Loop through fastq files
        fh, sp = self.__read_handle_create(fastq)
        for read in FastqGeneralIterator(fh):
            # Extract read name and number
            name, other = read[0].split(None, 1)
            number = other.split(':', 1)[0]
            # Add label to read name
            read = (
                '{}:{} {}'.format(name, number, other),
                read[1],
                read[2]
            )
            # Send data or break iteration
            try:
                self.__send((name, number, read), conn2)
            except StopIteration:
                break
        # Clean up
        if sp:
            sp.terminate()
        fh.close()
        conn2.close()
    
    def interleave_reads(
            self, outFastq
        ):
        ''' Function interleaves paired fastq files into a single fastq
        file.
        
        Args:
            label (bool)- Add ':1' label to read1 name and ':2' label to
                read2 name.
            check_pairs (bool)- Check read pairing prior to interleaving.
        
        Returns:
            count (int)- Number of paired reads processed
        
        '''
        # Check paired fastq files are present
        if not self.pair:
            raise ValueError('Paired fastq files required')
        # Check pairing if required
        self.check_names()
        # Loop through fastq files:
        for fastq in self.fastq_list:
            # Create pipe and process
            conn = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__interleave_read_process,
                args = (fastq, conn)
            )
            process.start()
            conn[1].close()
            # Store process data
            self.read_processes.append((process, conn[0]))
        # Extract data and count reads
        count = 0
        with writeFastq(outFastq, None, self.shell) as fastqOut:
            while True:
                try:
                    data = self.__read_processes_recv()
                except EOFError:
                    break
                count += 1
                # Check names and read number
                if (data[0][0] != data [1][0]
                    or data[0][1] != '1'
                    or data[1][1] != '2'):
                    self.__read_process_stop()
                    raise ValueError('Read {} name error'.format(count))
                # Write output fastq:
                fastqOut.write(data[0][2])
                fastqOut.write(data[1][2])
        # Stop read processes and return count
        self.__read_process_stop()
        return(count)
    
    def __interleave_trim_read_process(
            self, trim, fastq, conn
        ):
        ''' Function to trim fastq reads until after supplied sequence.
        
        Args:
            fastq (str)- Full path to the FASTQ file to read.
            trim (str)- Sequence after which sequence is trimmed.
            pend - End of multiprocessing pipe down which data will be sent.
        
        Sends:
            name (str)- Read name
            number (int)- Read number
        
        '''
        # Process pipes
        conn1, conn2 = conn
        conn1.close()
        # Extract numbe for string find adjustment
        adj = len(trim)
        # Loop through fastq files
        fh, sp = self.__read_handle_create(fastq)
        for read in FastqGeneralIterator(fh):
            # Extract read data
            header, sequence, quality = read
            name, description = header.split(None, 1)
            number = description.split(':', 1)[0]
            # Find trim sequence
            trimLoc = sequence.find(trim)
            if trimLoc != -1:
                trimLoc += adj
                sequence = sequence[:trimLoc]
                quality = quality[:trimLoc]
            # Send data or break iteration
            read = (
                '{}:{} {}'.format(name, number, description),
                sequence,
                quality
            )
            try:
                self.__send((name, number, trimLoc, read), conn2)
            except StopIteration:
                break
        # Clean up
        if sp:
            sp.terminate()
        fh.close()
        conn2.close()

    def interleave_trim_reads(
        self, trim, outFastq, minLength = 20
    ):
        ''' Function interleaves paired fastq files into a single fastq
        file.
        
        Args:
            label (bool)- Add ':1' label to read1 name and ':2' label to
                read2 name.
            check_pairs (bool)- Check read pairing prior to interleaving.
        
        Returns:
            count (int)- Number of paired reads processed
        
        '''
        # Check paired fastq files are present
        if not self.pair:
            raise ValueError('Paired fastq files required')
        # Loop through fastq files:
        for fastq in self.fastq_list:
            # Create pipe and process
            conn = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__interleave_trim_read_process,
                args = (trim, fastq, conn)
            )
            process.start()
            conn[1].close()
            # Store process data
            self.read_processes.append((process, conn[0]))
        # Extract data and count reads
        metrics = {'total':0, 'short':0, 'trim1':0, 'trim2':0}
        with writeFastq(outFastq, None, self.shell) as fastqOut:
            while True:
                try:
                    data = self.__read_processes_recv()
                except EOFError:
                    break
                metrics['total'] += 1
                # Check names and read number
                if data[0][0] != data [1][0]:
                    self.__read_process_stop()
                    raise ValueError('Read name mismatch')
                if (data[0][1] != '1'
                    or data[1][1] != '2'):
                    self.__read_process_stop()
                    raise ValueError('Read number error')
                # Count and skip short reads
                if (-1 <data[0][2] < minLength
                    or -1 < data[1][2] < minLength):
                    metrics['short'] += 1
                    continue
                # Write output fastq
                if data[0][2] != -1:
                    metrics['trim1'] += 1
                if data[1][2] != -1:
                    metrics['trim2'] += 1
                fastqOut.write(data[0][3])
                fastqOut.write(data[1][3])
        # Stop read processes and return count
        self.__read_process_stop()
        return(metrics)
    
class writeFastq(object):
    ''' An object that uses multiprocessing processes to parralelize the
    writing of FASTQ files. It assumes that the fastq files to write will
    be generated with the FastqGeneralIterator from the ??? package.
    Therefore read names will not have a leading '@' symbol which is
    therefore added by the write function.
    '''
    
    def __init__(self, fastq1, fastq2 = None, shell = True):
        ''' Function to initialise object. Function takes three arguments:
        
        1)  fastq1 - Full path to FASTQ file.
        2)  fastq2 - Full path to paired FASTQ file (optional).
        3)  shell - Boolean indicating whether to use shell gzip command
            to write gzipped output.
        
        '''
        # Store fastq files
        if not fastq2 is None:
            self.fastq_list = [fastq1, fastq2]
            self.pair = True
        else:
            self.fastq_list = [fastq1]
            self.pair = False
        ## Check output directories
        for fastq in self.fastq_list:
            outDir = os.path.dirname(fastq)
            if not os.path.isdir(outDir):
                raise IOError('Could not find output directory')
        # Store shell argument and process list
        self.shell = shell
        self.process_list = []
    
    def __write_process(self, fastq, conn):
        ''' Function to generate a process to write FASTQ files. FASTQ reads
        received from the multiprocessing pipe will be written to file. Receipt
        of None will cause the termination of the process. If self.shell is True
        then gzipped output files will be written using the gzip command in the
        shell. Funcion takes two arguments:
        
        1)  fastq - Full path to the FASTQ file to be created
        2)  pend - End of multiprocessing pipe from which reads will be
            extracted.
        '''
        # Process connections
        conn[1].close
        # Write gzip file using shell
        if self.shell and fastq.endswith('.gz'):
            # Create process
            command = 'gzip -c > %s' %(fastq)
            sp = subprocess.Popen(command, shell=True,
                stdin = subprocess.PIPE)
            fh = sp.stdin
        # Write file using python
        elif fastq.endswith('.gz'):
            sp = None
            fh = gzip.open(fastq, 'w')
        else:
            sp = None
            fh = open(fastq, 'w')
        # Extract reads from pipe and write to file
        while True:
            # Extract read
            read = conn[0].recv()
            # Break loop if read is none
            if read is None:
                break
            # Write read to file
            readString = '@{}\n{}\n+\n{}\n'.format(read[0], read[1], read[2])
            fh.write(readString)
        # Close files, pipes and subprocesses
        fh.close()
        if sp:
            sp.communicate()
        conn[0].close()
    
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
            conn = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = self.__write_process,
                args = (fastq, conn)
            )
            process.start()
            conn[0].close()
            # Store pipe end and processes
            self.process_list.append((process, conn[1]))
    
    def close(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. Thelist 
        self.process_list is then emptied.
        '''
        # Extract process and pipes
        for process, conn in self.process_list:
            # Add poisin pill, join process and close pipe
            try:
                conn.send(None)
            except IOError:
                pass
            process.join()
            conn.close()
        # Clear pipe and process list
        self.process_list = []
    
    def __pipe_send(self, read, conn):
        ''' Function to send a FASTQ read down a multiprocessing pipe. If
        the read does not correspond to the expected format all active
        processes are terminated and an IOError is raised. Function takes
        two arguments:
        
        1)  read - FASTQ read. This should consist of a tuple/list of
            three elements: read name, sequence, quality string.
        2)  pend - End of multiprocessing pipe down which read should be
            sent.
        '''
        # Add read1 to pipe
        if isinstance(read, (tuple,list)) and len(read) == 3:
            conn.send(read)
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
