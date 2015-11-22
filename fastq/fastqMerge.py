''' Functions to merge FASTQ files '''
import multiprocessing
import gzip
import re
from ngs_analysis.fastq import fastqExtract

def mergeLabelFastq(fastqIn1, fastqIn2, fastqOut, label1 = ':1',
    label2 = ':2'):
    ''' Function merges two FASTQ files into a single FASTQ file.
    Specified names are added to the end of the merged reads.
    
    1)  fastqIn1 - Python list containg read1 FASTQ files.
    2)  fastqIn2 - Python list containg read2 FASTQ files.
    3)  fastqOut - Output FASTQ file
    4)  label1 - Label to add to read1 file.
    5)  label2 - Label to add to read2 file.
    
    '''
    # Open output file
    if fastqOut.endswith('.gz'):
        outFile = gzip.open(fastqOut, 'w')
    else:
        outFile = open(fastqOut, 'w')
    # Extract reads from read1 and read2 files
    for f1, f2 in zip(fastqIn1, fastqIn2):
        # Create pipes
        pipe1Recv, pipe1Send = multiprocessing.Pipe(False)
        pipe2Recv, pipe2Send = multiprocessing.Pipe(False)
        # Create and start processes
        p1 = multiprocessing.Process(
            target = fastqExtract.read2Pipe,
            args = (f1, pipe1Send) 
        )
        p1.start()
        pipe1Send.close()
        p2 = multiprocessing.Process(
            target = fastqExtract.read2Pipe,
            args = (f2, pipe2Send) 
        )
        p2.start()
        pipe2Send.close()
        # Extract labelled reads and save to output
        while True:
            try:
                read1 = pipe1Recv.recv()
                read2 = pipe2Recv.recv()
            except EOFError:
                break
            # Extract read ID and check for equality
            read1ID = read1[0].split(' ' ,1)
            read2ID = read2[0].split(' ' ,1)
            if read1ID[0] != read2ID[0]:
                raise IOError('Files %s and %s contqain unmatched reads' %(
                    f1,f2))
            else:
                read1ID[0] += label1
                read2ID[0] += label2
                outFile.write('%s\n%s\n%s\n%s\n' %(
                    ' '.join(read1ID),
                    '\n'.join(read1[1:]),
                    ' '.join(read2ID),
                    '\n'.join(read2[1:])
                ))
    # Close pipes processes and files
    pipe1Recv.close()
    p1.join()
    pipe2Recv.close()
    p2.join()
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

