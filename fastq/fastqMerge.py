''' Functions to merge FASTQ files '''
import multiprocessing
import gzip
import re
from ngs_analysis.fastq import fastqExtract
from general_functions import writeFile

def mergeLabelPair(fastqIn1, fastqIn2, fastqOut, label1=':1', label2=':2'):
    ''' Function merges two FASTQ files into a single FASTQ file.
    Specified names are added to the end of the merged reads.
    
    1)  fastqIn1 - Read one FASTQ file(s). Either a string or a list of
        strings.
    2)  fastqIn2 - Read two FASTQ file(s). Either a string or a list of
        strings.
    3)  fastqOut - Output FASTQ file
    4)  label1 - Label to add to read1 file.
    5)  label2 - Label to add to read2 file.
    
    '''
    # Open input and output process
    outputProcess, outputPipe = writeFile.writeProcess(
        fileName = fastqOut, shell = True)
    input1Process, input1Pipe = fastqExtract.readToPipeProcess(fastqIn1)
    input2Process, input2Pipe = fastqExtract.readToPipeProcess(fastqIn2)
    # Extract labelled reads and save to output
    while True:
        try:
            read1 = input1Pipe.recv()
            read2 = input2Pipe.recv()
        except EOFError:
            break
        # Extract read ID and check for equality
        read1Header, read1Remainder = read1.split('\n' ,1)
        read1Header = read1Header.split(' ', 1)
        read2Header, read2Remainder = read2.split('\n' ,1)
        read2Header = read2Header.split(' ', 1)
        if read1Header[0] != read2Header[0]:
            raise IOError('Input FASTQ files contain unmatched reads')
        else:
            read1Header[0] += label1
            read2Header[0] += label2
            outputPipe.send('%s\n%s\n%s\n%s\n' %(
                ' '.join(read1Header),
                read1Remainder,
                ' '.join(read2Header),
                read2Remainder
            ))
    # Close pipes and processes
    input1Pipe.close()
    input2Pipe.close()
    input1Process.join()
    input2Process.join()
    outputPipe.close()
    outputProcess.join()

def mergeLabelTrimPair(fastqIn1, fastqIn2, trimSeq, fastqOut, minLength = 20,
    label1=':1', label2=':2'):
    ''' Function merges two paired FASTQ files into a single FASTQ file.
    FASTQ entries are trimmed to not extend beyong a supplied trim
    sequence. Any pair of reads for which one of the trimmed reads is
    shorter than the supplied minimum length is discarded. Specified
    labels are added to the end of the merged reads. Function takes seven
    arguments:
    
    1)  fastqIn1 - Read one FASTQ file(s). Either a string or a list of
        strings.
    2)  fastqIn2 - Read two FASTQ file(s). Either a string or a list of
        strings.
    3)  trimSeq - Sequence at which to terminate reads
    4)  fastqOut - Output FASTQ file.
    5)  minLength - Minimum length of reads to be output.
    6)  label1 - Label to add to read1 file.
    7)  label2 - Label to add to read2 file.
    
    Function returns a dictionary containing the following elements

    1)  total - Total number of read pairs.
    2)  short - Number of pairs with at least one read too short.
    3)  trim1 - Number of acceptable pairs with read1 trimmed.
    4)  trim2 - Number of acceptable pairs with read2 trimmed.
    '''
    # Open input and output process
    outputProcess, outputPipe = gzipFile.writeFromPipeProcess(
        fileName = fastqOut, shell = True)
    input1Process, input1Pipe = fastqExtract.readToPipeProcess(fastqIn1)
    input2Process, input2Pipe = fastqExtract.readToPipeProcess(fastqIn2)
    # Create output dictionary and key variables
    metrics = {'total' : 0, 'short': 0, 'trim1': 0, 'trim2' : 0}
    seqLength = len(trimSeq)
    # Extract labelled reads and save to output
    while True:
        try:
            read1 = input1Pipe.recv()
            read2 = input2Pipe.recv()
        except EOFError:
            break
        # Count total reads
        metrics['total'] += 1
        # Extract elements of read and identify trim sequence
        read1 = read1.split('\n')
        read2 = read2.split('\n')
        read1Loc = read1[1].find(trimSeq)
        read2Loc = read2[1].find(trimSeq)
        # Set trim length for read1 and count and skip if too short
        if read1Loc == -1:
            read1Trim = None
        else:
            read1Trim = read1Loc + seqLength
        if read1Trim and read1Trim < minLength:
            metrics['short'] += 1
            continue
        # Set trim length for read2 and count and skip if too short
        if read2Loc == -1:
            read2Trim = None
        else:
            read2Trim = read2Loc + seqLength
        if read2Trim and read2Trim < minLength:
            metrics['short'] +=1
            continue
        # Trim reads if required
        if read1Trim and read1Trim < len(read1[1]):
            read1[1] = read1[1][:read1Trim]
            read1[3] = read1[3][:read1Trim]
            metrics['trim1'] += 1
        if read2Trim and read2Trim < len(read2[1]):
            read2[1] = read2[1][:read2Trim]
            read2[3] = read2[3][:read2Trim]
            metrics['trim2'] += 1
        # Label reads
        if label1:
            read1Head = read1[0].split(' ',1)
            read1Head[0] += label1
            read1[0] = ' '.join(read1Head)
        if label2:
            read2Head = read2[0].split(' ',1)
            read2Head[0] += label2
            read2[0] = ' '.join(read2Head)
        # Save to file
        outputPipe.send('%s\n%s\n' %(
            '\n'.join(read1),
            '\n'.join(read2)
        ))
    # Close pipes and processes
    input1Pipe.close()
    input2Pipe.close()
    input1Process.join()
    input2Process.join()
    outputPipe.close()
    outputProcess.join()
    # Return metrics
    return(metrics)

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

