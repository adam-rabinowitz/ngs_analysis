import collections
import itertools
import os
import re

def parseModuleFile(modulefile):
    ''' Function to parse module file.
    
    Args:
        modulefile (str)- Path to tab-delimited file containing module data.
            The first and second columns are required and are the program and
            the program path. Additional columns should list the modules
            required for the program.
    
    Returns:
        pathDict (dict)- A dictionary where the key is the program and the
            value is the path
        moduleDict (dict)- A dictionary where the key is the program and the
            value is a list of required modules.
    '''
    # Create output variables
    pathDict = {}
    moduleDict = {}
    # Import each line as a list
    with open(modulefile) as infile:
        for line in infile:
            linedata = line.strip().split('\t')
            # Extract and store data
            program = linedata[0]
            path = linedata[1]
            modules = linedata[2:]
            pathDict[program] = path
            moduleDict[program] = modules
    # Return data
    return(pathDict, moduleDict)

def parseSampleFile(samplefile):
    ''' Function to parse sample file.
    
    Args:
        samplefile (str)- Path to tab-delimited file containing sample data.
            The first column is the sample name which will be used as a prefix
            for all outut files. The second column should be the prefix for
            the identification of FASTQ files. Additional columns should list
            directories in which to search for FASTQ files.
    
    Returns:
        sampleDict (dict)- A collections ordered dictionary where the key
            is the sample name and the value in a tuple where the first
            element is the prefix and the second element is a list of
            directories.
    
    '''
    # Create output variable
    sampleDict = collections.OrderedDict()
    prefixList = []
    # Import each line of the file as list
    with open(samplefile) as infile:
        for line in infile:
            linedata = line.strip().split('\t')
            # Extract and store data
            name = linedata[0]
            prefix = linedata[1]
            indirs = linedata[2:]
            if len(indirs) < 1:
                raise IOError('No input directores for {}'.format(name))
            sampleDict[name] = (prefix, indirs)
            prefixList.append(prefix)
    # Check prefixes will identify unique files
    for p1, p2 in itertools.permutations(prefixList, 2):
        if p1.startswith(p2):
            raise IOError("prefices '{}' and '{}' overlap".format(p1, p2))
    # Return output
    return(sampleDict)

def parseParameterFile(paramfile):
    ''' Function to parse parameter file
    
    Args:
        paramfile (str)- Path to tabdelimited paramter file.
    
    Returns:
        paramDict (dict)- An dictionary of all parameters for analysis.
    
    '''
    # Create and populate parameter file
    paramDict = {}
    with open(paramfile) as infile:
        for line in infile:
            # Skip comment lines
            if line.startswith('#'):
                continue
            # Extract data and try type conversion
            param, value = line.strip().split('\t')
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
            # Store data
            paramDict[param] = value
    # Return data
    return(paramDict)

def parseIndexFile(indexfile):
    ''' Function to parse index file
    
    Args:
        paramfile (str)- Path to tabdelimited paramter file.
    
    Returns:
        paramDict (dict)- An dictionary of all parameters for analysis.
    
    '''
    # Create and populate index dictionary
    indexDict = {}
    with open(indexfile) as infile:
        for line in infile:
            # Extract data and try type conversion
            param, value = line.strip().split('\t')
            indexDict[param] = value
    # Return data
    return(indexDict)

def findFastq(prefix, dirList):
    ''' A function to identify FASTQ files from directories
    using a supplied filename prefix
    
    Args:
        prefix (str)- Prefix of the FASTQ files to be found.
        dirList (list)- A list of directories to search.
        
    Returns:
        read1 (list)- A list of read1 FASTQ files
        read2 (list)- A list of read2 FASTQ files
    '''
    # Create variables to store results
    read1 = []
    read2 = []
    # Create regular expression to find files
    prefix = re.escape(prefix)
    read1Pattern = re.compile(prefix + '.*?R1(_\\d{3}){0,1}\\.fastq.gz$')
    # Loop through directories to find fastq files
    for directory in dirList:
        # Loop through file names and find read1 files
        filenames = os.listdir(directory)
        for f in filenames:
            if re.match(read1Pattern, f):
                read1.append(os.path.join(directory, f))
                # Find and store matching read2 files
                read2File, nsub = re.subn(
                    'R1(?=(_\\d{3}){0,1}\\.fastq.gz$)', 'R2', f)
                if nsub != 1:
                    raise IOError('Could not generate read2 filename'\
                        ' for %s' %(f))
                if read2File in filenames:
                    read2.append(os.path.join(directory, read2File))
    # Check output files and return
    if len(read1) == 0:
        raise IOError('{}: No FASTQ files found'.format(prefix))
    if len(read2) and len(read1) != len(read2):
        raise IOError('{}: Mixed single- and paired-end'.format(prefix))
    return(read1, read2)

def createOutFiles(outdir, sample):
    ''' Function to create output files for analysis
    
    Args:
        outdir (str)- Path to output directory
        sample (str)- Sample name
    
    Returns
        outDict (dict)- Dictionary of output files.
    
    '''
    # Create variable to store files
    outfiles = {}
    # Create output directories and output prefix
    sampledir = os.path.join(outdir, sample)
    if not os.path.isdir(sampledir):
        os.mkdir(sampledir)
    outprefix = os.path.join(sampledir, sample) + '.'
    # Store directories, prefixes and job file
    outfiles['prefix'] = outprefix
    outfiles['outdir'] = sampledir
    outfiles['slurm'] = outprefix + 'slurm'
    # Create file names for processing FASTQ files
    outfiles['cat1'] = outprefix + 'R1.fastq.gz'
    outfiles['cat2'] = outprefix + 'R2.fastq.gz'
    outfiles['trim1'] = outprefix + 'trim.R1.fastq.gz'
    outfiles['trim2'] = outprefix + 'trim.R2.fastq.gz'
    outfiles['fastqclog'] = outprefix + 'fastqc.log'
    outfiles['trimlog'] = outprefix + 'cutadapt.metrics'
    # Create file names for processing BAM files
    outfiles['starbam'] = outprefix + 'Aligned.out.bam'
    outfiles['starlog'] = outprefix + 'star.log'
    outfiles['sortbam'] = outprefix + 'sort.bam'
    outfiles['sortlog'] = outprefix + 'sort.log'
    outfiles['mdupbam'] = outprefix + 'mdup.bam'
    outfiles['mduplog1'] = outprefix + 'mdup.metrics'
    outfiles['mduplog2'] = outprefix + 'mdup.log'
    # Create output files for htseq
    outfiles['htseqlog'] = outprefix + 'htseq.log'
    outfiles['genecounts'] = outprefix + 'gene_counts.txt'
    # Create file names for QC of BAM files
    outfiles['metrlog1'] = outprefix + 'collectrna.metrics'
    outfiles['metrlog2'] = outprefix + 'collectrna.log'
    outfiles['alsumlog1'] = outprefix + 'alignsum.metrics'
    outfiles['alsumlog2'] = outprefix + 'alignsum.log'
    # Return data
    return(outfiles)
    
def fastQC(inFile, outDir, path):
    ''' This function performs a FastQC analysis on a fastq file and
    then organises the output data. Function is built for version 0.11.2
    of FastQC. Function takes three arguments:
    
    1)  inFile - Input FASTQ file.
    2)  outDir - Output directory.
    3)  path - Path to FastQC; Default = 'fastqc'.
    
    '''
    # Extract sample name
    name = re.search('([^/]+)\\.fastq(?:\\.gz){0,1}$',inFile).group(1)
    # Create FastQC command and return it
    fastqcCommand = '%s --extract -q -o %s %s && rm %s %s' %(
        path,
        outDir,
        inFile,
        os.path.join(outDir, name + '_fastqc.html'),
        os.path.join(outDir, name + '_fastqc.zip')
    )
    # Execute or return command
    return(fastqcCommand)

def cutadapt(
        read1In, read1Out, read2In, read2Out, quality, adapter, length, path,
        overlap, error
    ):
    ''' A function to create cutadapt command
    
    Args:
        read1In (str)- Path to read1 input file.
        read1Out (str)- Path to read2 output file.
        read2In (str)- Path to read2 input file.
        read2Out (str)- Path to read2 output file.
        quality (int)- Base quality score to use for trimming.
        adapter (str)- Adapter to use for trimming.
        length (int)- Minimum length of trimmed reads.
        path (str)- Path for cutadapt program.
    
    '''
    # Check arguments
    if not read2In is None and read2Out is None:
        raise IOError('Output file must be supplied for 2nd read')
    if not isinstance(length, int):
        raise TypeError('length must be integer')
    if length < 25:
        raise ValueError('length must be >=25')
    if not isinstance(overlap, int):
        raise TypeError('overlap must be integer')
    if not 1 <= overlap <= len(adapter):
        raise ValueError('overlap must be >=1 and <= adapter length')
    if not isinstance(error, (int, float)):
        raise TypeError('error must be integer or float')
    if not 0 <= error < 1:
        raise ValueError('error must be >=0 and <1')
    # Create single end argument
    adapterList = adapter.split(',')
    command = [path]
    if read2In is None:
        for a in adapterList:
            command.extend(['-a', a])
        command.extend([
            '-o', read1Out, '-e', error, '-q', quality, '-m', length, '-O',
            overlap, read1In])
    else:
        for a in adapterList:
            command.extend(['-a', a, '-A', a])
        command.extend([
            '-o', read1Out, '-p', read2Out, '-e', error, '-q', quality, '-m',
            length, '-O', overlap, read1In, read2In])
    # Join and return command
    command = ' '.join(map(str, command))
    return command

def starAlign(
        indexDir, outPrefix, read1, read2, threads, path, rg=1,
        pl='uknown', lb='unknown', sm='uknown'
    ):
    # Create output command
    command = [path, '--runThreadN', threads, '--genomeDir', indexDir,
        '--outFileNamePrefix', outPrefix, '--outSAMtype', 'BAM', 'Unsorted',
        '--outSAMunmapped', 'Within', '--readFilesIn', read1]
    if read2:
        command.append(read2)
    # Append read file command
    if read1.endswith('.gz'):
        if read2.endswith('.gz'):
            command.extend(['--readFilesCommand', 'zcat'])
        else:
            raise ValueError('mixture of compressed and uncompressed files')
    # Add read group information
    if rg:
        command.extend(['--outSAMattrRGline', 'ID:{}'.format(rg)])
        if pl:
            command.append('PL:{}'.format(pl))
        if lb:
            command.append('LB:{}'.format(lb))
        if sm:
            command.append('SM:{}'.format(sm))
    # Concatenate commadn and return
    command = ' '.join(map(str, command))
    return(command)

def bamsort(
        inFile, outFile, threads, memory, path
    ):
    ''' Function to create sort BAM commabd using samtools.
    
    Args:
        inFile (str)- Path to input file.
        outFile (str)- Path to outfile.
        threads (int)- Number of threads to use in sort.
        memory (int)- Memory, in gigabytes, to use in each thread.
        path (str)- Path to samtools executable.
    
    Returns:
        sortCommand (str)- Output command
    
    '''
    # Check input file
    if not inFile.endswith('.bam'):
        raise TypeError('Input file suffix must be .bam')
    # Check output file
    if not outFile.endswith('.bam'):
        raise TypeError('Output file suffix must be .bam')
    # Process memory argument
    memory = str(memory) + 'G'
    # Generate sort command
    sortCommand = [path, 'sort', '-m', memory, '-@', str(threads),
        '-o', outFile, '-T', outFile[:-4], '-O', 'BAM', inFile]
    sortCommand = filter(None, sortCommand)
    sortCommand = ' '.join(sortCommand)
    # Delete input and index output
    sortCommand += ' && {} index {}'.format(path, outFile)
    sortCommand += ' && rm {}'.format(inFile)
    # Return command
    return(sortCommand)

def markDuplicates(
        inBam, outBam, logFile, picardPath, memory
    ):
    ''' Function to mark duplicates using the picard toolkit.
    
    Args:
        inBam (str)- Full path to input BAM file.
        outBam (str)- Full path to output BAM file.
        logFile (str)- Full path to output log file.
        picardPath (str)- Path to picard jar file.
        memory (int)- Amount of memory in java heap in gigabytes.
    
    Returns:
        command (str)- Mark duplicates command
    
    '''
    # Create command
    command = [
        'java', '-jar', '-Xmx{}g'.format(memory), picardPath, 'MarkDuplicates',
        'I=' + inBam, 'O=' + outBam, 'M=' + logFile, 'ASSUME_SORTED=true',
        'CREATE_INDEX=true', 'REMOVE_DUPLICATES=false'
    ]
    # Merge command, add deletion and return
    command = ' '.join(command)
    command += ' && rm {}*'.format(inBam[:-1])
    return(command)

def rnaseqMetric(
        bam, output, refflat, strand, rrna, path, memory
    ):
    ''' Function to generate command for picard CollectRNASeqMetrics
    
    Args:
        bam (str)- Path to input BAM file.
        output (str)- Path to output file.
        refflat (str)- Path to reflat file.
        strand (str)- Strand: should be one none|forward|reverse.
        
     Returns:
        command (str)- CollectRnaSeqMetrics command.
    
    '''
    # Check strand argument
    if strand == 'none':
        strandArg = 'STRAND=NONE'
    elif strand == 'forward':
        strandArg = 'STRAND=FIRST_READ_TRANSCRIPTION_STRAND'
    elif strand == 'reverse':
        strandArg = 'STRAND=SECOND_READ_TRANSCRIPTION_STRAND'
    else:
        raise ValueError('strans must be one of none|forward|reverse')
    # Build command
    command = [
        'java', '-jar', '-Xmx{}g'.format(memory), path, 'CollectRnaSeqMetrics',
        'I=' + bam, 'O=' + output, 'REF_FLAT=' + refflat, strandArg,
        'RIBOSOMAL_INTERVALS=' + rrna
    ]
    # Join and return command
    command = ' '.join(command)
    return(command)

def alignMetrics(
        bam, output, fasta, path, memory
    ):
    ''' Function to generate command for picard CollectAlignmentSummeryMetrics
    
    Args:
        bam (str)- Path to input BAM file.
        output (str)- Path to output file.
        fasta (str)- Path to FASTA file.
        path (str)- Path to picard executable file.
        memory (int)- Initial heap size in gigabytes.
    
    Returns:
        command (str)- CollectAlignmentSummaryMetrics command.
    
    '''
    # Create command
    command = [
        'java', '-jar', '-Xmx{}g'.format(memory), path,
        'CollectAlignmentSummaryMetrics', 'R=' + fasta, 'I=' + bam,
        'O=' + output
    ]
    # Join and return command
    command = ' '.join(command)
    return(command)

def htseq(
        bam, gtf, path, feature='exon', attrid='gene_id', mode='union',
        stranded='reverse', mapq=10
    ):
    # Check arguments
    if not mode in ('union', 'intersection-strict', 'intersection-nonempty'):
        raise ValueError('unrecognised mode')
    if not stranded in ('yes', 'no', 'reverse'):
        raise ValueError('unrecognised stranded argument')
    if not isinstance(mapq, int):
        raise TypeError('mapq not an integer')
    if mapq < 0:
        raise ValueError('mapq is negative')
    # Create command
    command = [path, '-f', 'bam', '-r', 'pos', '-s', stranded, '-t', feature,
        '-i', attrid, '-m', mode, '-a', mapq, bam, gtf]
    # Join and return command
    command = ' '.join(map(str, command))
    return(command)
