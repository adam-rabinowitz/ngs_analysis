################################################################################
## cutadaptTrimPaired
################################################################################
# Define function
def cutadaptTrimPaired(read1In, read2In, read1Out, read2Out, quality = 20,
    adapter = 'AGATCGGAAGAGC', length = 25, path = 'cutadapt'
):
    ''' This function peforms quality trimming on paired end FASTQ files
    using the cutadapt package. Function is built for version 1.7.1 of
    cutadapt which uses a two stage process to trim FASTQ file. The
    default adapter sequence is common to both ends of the standard
    Illumina adapter. Output file names ending in '.gz' will cause the
    output to be gzipped compressed. Function takes  8 arguments:
    
    1)  Read1 input fastq file or list of FASTQ files.
    2)  Read2 input fastq file or list of FASTQ files.
    3)  Read1 output fastq file.
    4)  Read2 output fastq file.
    5)  Minimum base quality; Default = 20.
    6)  Adapter sequence; Default = 'AGATCGGAAGAGC'.
    7)  Minimum read length; Default = 25. 
    8)  Path to cutadapt; Default = 'cutadapt'.
    
    '''
    # Check name of output files
    if (not read1Out.endswith('.fastq.gz') or
        not read2Out.endswith('.fastq.gz')):
        raise IOError('Output must be gzipped FASTQ files')
    # Create temporary FASTQ files
    concatRead1 = read1Out[:-9] + '.concat.fastq.gz'
    concatRead2 = read2Out[:-9] + '.concat.fastq.gz'
    tempRead1 = read1Out[:-9] + '.trimtmp.fastq.gz'
    tempRead2 = read2Out[:-9] + '.trimtmp.fastq.gz'
    # Process read1In list input
    if isinstance(read1In, list):
        # Check suffixes:
        for read in read1In:
            if not read.endswith('.fastq.gz'):
                raise IOError('Input must be gzipped FASTQ files')
        # Create concatenation command
        if len(read1In) > 1:
            concat1Command = 'zcat %s | gzip > %s' %(
                ' '.join(read1In), concatRead1)
            trimFile1 = concatRead1
        else:
            concat1Command = ''
            trimFile1 = read1In[0]
    # Process read1In string input
    elif isinstance(read1In, str):
        if not read1In.endswith('.fastq.gz'):
            raise IOError('Input must be gzipped FASTQ files')
        concat1Command = ''
        trimFile1 = read1In
    # Process read2In list input
    if isinstance(read2In, list):
        # Check suffixes:
        for read in read2In:
            if not read.endswith('.fastq.gz'):
                raise IOError('Input must be gzipped FASTQ files')
        # Create concatenation command
        if len(read2In) > 1:
            concat2Command = 'zcat %s | gzip > %s' %(
                ' '.join(read2In), concatRead2)
            trimFile2 = concatRead2
        else:
            concat2Command = ''
            trimFile2 = read2In[0]
    # Process read1In string input
    elif isinstance(read2In, str):
        if not read2In.endswith('.fastq.gz'):
            raise IOError('Input must be gzipped FASTQ files')
        concat2Command = ''
        trimFile2 = read2In
    # Generate trimming commands
    trim1Command = [path, '-q', str(quality), '-a', adapter, '-e',
        '0.1', '--minimum-length', str(length), '-O', '1', '-o', tempRead1,
        '-p', tempRead2, trimFile1, trimFile2]
    trim2Command = [path, '-q', str(quality), '-a', adapter, '-e',
        '0.1', '--minimum-length', str(length), '-O', '1', '-o', read2Out,
        '-p', read1Out, tempRead2, tempRead1]
    # Generate command to remove temporary files
    removeCommand = 'rm %s %s' %(tempRead1, tempRead2)
    if concat1Command:
        removeCommand += ' %s' %(concatRead1)
    if concat2Command:
        removeCommand += ' %s' %(concatRead2)
    # Join and return command
    jointCommand = [concat1Command, concat2Command, ' '.join(trim1Command),
        ' '.join(trim2Command), removeCommand]
    jointCommand = filter(None, jointCommand)
    jointCommand = " && ".join(jointCommand)
    return(jointCommand)



################################################################################
## cutadaptTrim
################################################################################
# Define function
def cutadaptTrim(
        readIn, readOut, quality = 20, adapter = 'AGATCGGAAGAGC', length = 25,
        path = 'cutadapt'
    ):
    ''' This function peforms quality trimming on a FASTQ file using the
    cutadapt package. Function is built for version 1.7.1 of cutadapt. The
    default adapter sequence is common to both ends of the standard
    Illumina adapter. Output file names ending in '.gz' will cause the
    output to be gzipped compressed. Function takes 6 arguments:
    
    1)  readIn - Input fastq file.
    2)  readOut - Output fastq file.
    3)  quality - Minimum base quality.
    4)  adapter - Adapter sequence to remove from 3' end of reads.
    5)  length - Minimum read length post trimming. 
    6)  path - Path to cutadapte executable.
    
    '''
    # Check name of output files
    if not readOut.endswith('.fastq.gz'):
        raise IOError('Output must be gzipped FASTQ file')
    # Create temporary FASTQ file
    concatRead = readOut[:-9] + '.concat.fastq.gz'
    # Process readIn list input
    if isinstance(readIn, list):
        # Check suffixes:
        for read in readIn:
            if not read.endswith('.fastq.gz'):
                raise IOError('Input must be gzipped FASTQ files')
        # Create concatenation command
        if len(readIn) > 1:
            concatCommand = 'zcat %s | gzip > %s' %(
                ' '.join(readIn), concatRead)
            trimFile = concatRead
        else:
            concatCommand = ''
            trimFile = readIn[0]
    # Process readIn string input
    elif isinstance(readIn, str):
        if not readIn.endswith('.fastq.gz'):
            raise IOError('Input must be gzipped FASTQ files')
        concatCommand = ''
        trimFile = readIn
    # Generate trimming commands
    trimCommand = [path, '-q', str(quality), '-a', adapter, '-e',
        '0.1', '--minimum-length', str(length), '-O', '1', '-o', readOut,
        trimFile]
    # Generate command to remove temporary files
    if concatCommand:
        removeCommand = 'rm %s' %(concatRead)
    else:
        removeCommand = ''
    # Join and return command
    jointCommand = [concatCommand, ' '.join(trimCommand), removeCommand]
    jointCommand = filter(None, jointCommand)
    jointCommand = " && ".join(jointCommand)
    return(jointCommand)
