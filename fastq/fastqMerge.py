''' Functions to merge FASTQ files '''
import multiprocessing
import gzip
from Bio import SeqIO
from github.general_functions import iohandle

def labelFastq(fastq, label, output):
    ''' Function to open FASTQ file, rename entries and return reads.
    Fnction takes three arguments:
    
    1)  fastq - Full path to fastq file.
    2)  label - Label to add to FASTQ name.
    3)  output - Output object handled by the iohandle.hadleout function.
    
    '''
    # Create object to handle output
    outObject = iohandle.handleout(output)
    # Open FASTQ file
    if fastq.endswith('.gz'):
        fastqHandle = gzip.open(fastq, 'r')
    else:
        fastqHandle = open(fastq, 'r')
    # Loop through reads, add label and add to output
    readData = SeqIO.parse(fastqHandle, 'fastq')
    for read in readData:
        read.id = read.id + label
        outObject.add(read.format('fastq'))
    # Close and return output object
    return(outObject.close())

def mergeLabelFastq(read1, read2, outFastq):
    ''' Function merges two FASTQ files into a single FASTQ file. A ':1'
    label is added to the end of read1 names and a ':2' label is added to
    the end of read2 names. Function takes three arguments:
    
    1)  read1 - Python list containg full path to read1 input files.
    2)  read2 - Python list containg full path to read2 input files.
    3)  outFastq - Full path to output file
    
    '''
    # Open output file
    if outFastq.endswith('.gz'):
        outFile = gzip.open(outFastq, 'w')
    else:
        outFile = open(outFastq, 'w')
    # Loop though read1 and read2 files
    for r1, r2 in zip(read1,read2):
        # Create pipes
        pipe1Recv, pipe1Send = multiprocessing.Pipe(False)
        pipe2Recv, pipe2Send = multiprocessing.Pipe(False)
        # Create and start processes
        p1 = multiprocessing.Process(
            target = labelFastq,
            args = (r1, ':1', pipe1Send) 
        )
        p1.start()
        pipe1Send.close()
        p2 = multiprocessing.Process(
            target = labelFastq,
            args = (r2, ':2', pipe2Send) 
        )
        p2.start()
        pipe2Send.close()
        # Extract labelled reads and save to output
        while True:
            try:
                outFile.write(pipe1Recv.recv() + pipe2Recv.recv())
            except EOFError:
                pipe1Recv.close()
                p1.join()
                pipe2Recv.close()
                p2.join()
                break
    # Close output file
    outFile.close()

def concatFastq(inPrefix, inDirList, outPrefix, pair=True):
    ''' Function to generate command to concatenate fastq files '''
    # Extract reads
    reads = findIlluminaFastq(prefix = inPrefix, dirList = inDirList,
        pair = pair)
    commands = []
    # Sequentially process read lists
    for readNo, readList in enumerate(reads):
        # Count gzip files
        gzipCount = sum([x.endswith('.gz') for x in readList])
        # Check that there are multiple files to be concatendated
        if len(readList) < 2:
            print "Prefix %s, %s, has %s files and won't be processed" %(
                inPrefix,
                'R' + str(readNo + 1),
                len(readList)
            )
            command = ''
        # Generate commands to perform concatenation
        elif gzipCount == 0:
            command = 'cat %s | gzip > %s' %(
                ' '.join(readList),
                outPrefix + '_R' + str(readNo + 1) + '.fastq.gz'
            )
        elif gzipCount == len(readList):
            command = 'zcat %s | gzip > %s' %(
                ' '.join(readList),
                outPrefix + '_R' + str(readNo + 1) + '.fastq.gz'
            )
        else:
            raise TypeError('Mix of unzipped and gzipped files')
        commands.append(command)
     # Return commands
    return(commands)

