from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

def read2Pipe(fastqFile, pipe):
    ''' Function parses a FASTQ file using Bio.SeqIO and creates a list
    object where each element is a line of the fastq file. The list elements
    are sent down the pipe.
            
    '''
    # Create command to open file
    if fastqFile.endswith('.gz'):
        openFile = gzip.open
    else:
        openFile = open
    # Create output
    for title, seq, qual in FastqGeneralIterator(openFile(fastqFile)):
        pipe.send([title,seq,'+',qual])
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
