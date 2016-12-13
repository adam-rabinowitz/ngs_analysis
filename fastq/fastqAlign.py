import os
import sys
from ngs_python.bam import samtools

def bowtie2Align(
        index, outFile, read1, read2 = None, bowtie2Path = 'bowtie2',
        threads = 1, readGroup = 1, sampleName = None, libraryID = None,
        platform = None, discordant = False, mixed = False, upto = None,
        maxInsert = None, check = True, samtoolsPath = 'samtools',
        memory = '2', nameSort = False
    ):
    ''' Function to generate command to peform Bowtie2 Alignment of
    paired FASTQ files. Function takes 9 arguments:
    
    1)  index - Suffix of Bowtie2 index.
    2)  outSam - Name of output SAM file.
    3)  read1 - Read1 FASTQ file.
    4)  read2 - Read2 FASTQ file.
    5)  path - Path to Bowtie2 executable.
    6)  threads - Number of thread to use.
    7)  discordant - Boolean; whether to output discordant pairs.
    8)  mixed - Boolean; whether to output mixed pairs.
    9)  upto - Number of reads to align
    10) check - Boolean; whether to check for index entensions.
    
    '''
    # Check for index extensions
    if not isinstance(check, bool):
        raise TypeError('check argument must be bool')
    if check:
        suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
            '.rev.2.bt2']
        for s in suffixes:
            if not os.path.isfile(index + s):
                raise IOError('Index file %s no found' %(index + s))
    # Check and process discordant
    if not isinstance(discordant, bool):
        raise TypeError('discordant argument must be bool')
    if discordant:
        discordant = ''
    else:
        discordant = '--no-discordant'
    # Check mixed
    if not isinstance(mixed, bool):
        raise TypeError('mixed argument must be bool')
    if mixed:
        mixed = ''
    else:
        mixed = '--no-mixed'
    # Check upto argument
    if not upto is None:
        if not isinstance(upto, int):
            raise TypeError('upto argument must be integer')
        if upto < 1:
            raise ValueError('upto argument must be >= 1')
    # Check maximum insert argument
    if not maxInsert is None:
        if not isinstance(maxInsert, int):
            raise TypeError('maxInsert argument must be integer')
        if maxInsert < 1:
            raise ValueError('maxInsert argument must be >= 1')
    # Check outut file name and generate intermediate file names
    if outFile.endswith('.sam'):
        outSam = outFile
        outBam = ''
    elif outFile.endswith('.bam'):
        outBam = outFile
        outSam = outFile[:-4] + '.sam'
    else:
        raise ValueError("'outFile' argument must end '.sam' or '.bam'")
    # Join multiple fastq files
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
    # Create initial command
    bowtie2Command = [bowtie2Path, '--phred33', '--very-sensitive', mixed,
        discordant, '-p', str(threads), '-x', index, '-S', outSam]
    # Extend command depending on if read2 is applied
    if read2:
        bowtie2Command.extend(['-1', read1, '-2', read2])
    else:
        bowtie2Command.extend(['-U', read1])
    # Supplement additional commands
    if upto:
        bowtie2Command.extend(['-u', str(upto)])
    if maxInsert:
        bowtie2Command.extend(['-X', str(maxInsert)])
    # Add read group data
    if readGroup:
        # Create read group list
        rgList = ['--rg-id', str(readGroup)]
        if sampleName:
            rgList.extend(['--rg', 'SM:' + str(sampleName)])
        if libraryID:
            rgList.extend(['--rg', 'LB:' + str(libraryID)])
        if platform:
            rgList.extend(['--rg', 'PL:' + str(platform)])
        # Add list to command
        bowtie2Command.extend(rgList)
    # Concatenate bowtie2Command command
    bowtie2Command = filter(None, bowtie2Command)
    bowtie2Command = ' '.join(bowtie2Command)
    # Supplement BWA command with sort command
    if outBam:
        sortCommand = samtools.sort(inFile = outSam, outFile = outBam,
            name = nameSort, memory = memory, delete = True,
            path = samtoolsPath, threads = threads)
        completeCommand = bowtie2Command + ' && ' + sortCommand
    else:
        completeCommand = bowtie2Command
    # Return complete command
    return(completeCommand)

def bwaMemAlign(
        index, outFile, read1, read2 = None, bwaPath = 'bwa', threads = 1,
        readGroup = '1', sampleName = None, libraryID = None, platform = None,
        markSecondary = True, check = True, samtoolsPath = 'samtools',
        memory = '2', nameSort = False
    ):
    ''' Function to generate command to perform BWA mem alignment of single
    end or paired end FASTQ files. If the supplied output file name ends with
    '.bam' then a sorted BAM file will be generated else if the file names
    ends with '.sam' a sam file is returned. Function takes the following 14
    arguments:
    
    Args:
        index (str)- Full path BWA index prefix.
        outFile (str)- Full path to output sam or bam file.
        read1 (str)- Read1 FASTQ file.
        read2 (str)- Read2 FASTQ file.
        bwaPath (str)- BWA exectuable.
        threads (int)- Number of threads to use.
        readGroup (str)- Read group to be used in SAM/BAM.
        sampleName (str)- Name of sample to be used in header.
        libraryID (str)- Library ID to be used in SAM/BAM.
        platform (str)- Platform to be used in SAM/BAM.
        markSecondary (bool)- Mark secondary alignments.
        check (bool)- Check for index extensions and output directory.
        samtoolsPath (str)- Samtools executable.
        memory (int)- Gigabytes of memory to use in generating BAM file.
        nameSort (bool)- Generate a name sorted BAM file.
    
    Returns:
        bwaCommand (str)- Command to perform BWA alignment.
    
    Raises:
        IOError - If index suffixes or output directory are absent.
        TypeError - If arguments are of the wrong type.
        ValueError - If arguments have an unexpected value.
    
    '''
    # Check index extensions and output directory, if required
    if not isinstance(check, bool):
        raise TypeError('check argument must be bool')
    if check:
        suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        for s in suffixes:
            if not os.path.isfile(index + s):
                raise IOError('Genome index file %s no found' %(index + s))
        if not os.path.isdir(os.path.dirname(outFile)):
            raise IOError('Could not find output directory {}'.format(
                os.path.dirname(outFile)))
    # Check outut file name and generate intermediate file names
    if outFile.endswith('.sam'):
        outSam = outFile
        outBam = ''
    elif outFile.endswith('.bam'):
        outBam = outFile
        outSam = outFile[:-4] + '.sam'
    else:
        raise ValueError('outFile argument must end .sam or .bam')
    # Process secondary command
    if not isinstance(markSecondary, bool):
        raise TypeError('markSecondary argument must be bool')
    if markSecondary:
        markSecondary = '-M'
    else:
        markSecondary = ''
    # Create command
    bwaCommand = [bwaPath, 'mem', markSecondary ,'-t', str(threads),
        index, read1, read2]
    # Remove missing elements from coomand
    bwaCommand = filter(None, bwaCommand)
    # Add read group data
    if readGroup:
        # Create read group string
        rgString = "'@RG\\tID:" + str(readGroup)
        if sampleName:
            rgString += '\\tSM:' + str(sampleName)
        if libraryID:
            rgString += '\\tLB:' + str(libraryID)
        if platform:
            rgString += '\\tPL:' + str(platform)
        rgString += "'"
        # Add string to command
        bwaCommand.insert(2,rgString)
        bwaCommand.insert(2,'-R')
    # Complete BWA command
    bwaCommand = '%s > %s' %(' '.join(bwaCommand), outSam)
    # Supplement BWA command with sort command, if required, and return
    if outBam:
        sortCommand = samtools.sort(inFile = outSam, outFile = outBam,
            name = nameSort, memory = memory, delete = True,
            path = samtoolsPath, threads = threads)
        bwaCommand = bwaCommand + ' && ' + sortCommand
    return(bwaCommand)

def rsemBowtie2Align(
        index, outPrefix, read1, read2 = None,
        rsemPath = 'rsem-calculate-expression', bowtie2Path = '', threads = 1,
        forProb = 0.5, genomeBam = True, estimateRspd = True, check = True
    ):
    ''' Function generates and returns a command to align paired FASTQ
    files to a Bowtie2 indexed RSEM refrence transcriptome.
    
    Args:
        index (str)- Full path to transcriptome index prefix.
        outPrefix (str)- Full path for output files.
        read1 (str/list)- Full path to read1 FASTQ file or list of paths.
        read2 (str/list)- Full path to read2 FASTQ file or list of paths.
        rsemPath (str)- Full path to rsem-calculate-expression executable.
        bowtie2Path (str)- Full path to directory containing bowtie2
            executables.
        threads (int)- Number of threads to use in alignment
        forProb (int/float)- Probabiliy that the first read is derived from
            the forward strand of the transcript.
        genomeBam (bool)- Generate a genome BAM file.
        esitmateRspd (bool)- Estimate read start position distribution.
        check (bool)- Check if index files are present.
    
    Returns:
        rsemCommand (str)- Command to perform RSEM alignment
    
    Raises:
        TypeError - If arguments are of the wrong type.
        IOError - If index files or output directory are absent.
        ValueError - If arguments are of unexpected values.
    
    '''
    # Check for correct suffixes and output directory
    if not isinstance(check, bool):
        raise TypeError('check argument must be bool')
    if check:
        suffix = [
            '.grp', '.ti', '.transcripts.fa', '.seq', '.idx.fa', '.1.bt2',
            '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'
        ]
        for s in suffix:
            if not os.path.isfile(index + s):
                raise IOError('Could not find index file {}'.format(
                    index + s))
        if not os.path.isdir(os.path.dirname(outPrefix)):
            raise IOError('Could not find output directory {}'.format(
                os.path.dirname(outPrefix)))
    # Process genomeBAM argument
    if not isinstance(genomeBam, bool):
        raise TypeError('genomeBam argument must be bool')
    if genomeBam:
        genomeBam = '--output-genome-bam'
    else:
        genomeBam = ''
    # Process estimateRspd argument
    if not isinstance(estimateRspd, bool):
        raise TypeError('estimateRspd argument must be bool')
    if estimateRspd:
        estimateRspd = '--estimate-rspd'
    else:
        estimateRspd = ''
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if not isinstance(read1, str):
        raise TypeError('read1 must be string or list of strings')
    if not read2 is None:
        if isinstance(read2, list):
            read2 = ','.join(read2)
        if not isinstance(read2, str):
            raise TypeError('read2 must be string or list of strings')
    # Build command
    rsemCommand = [rsemPath, '-p', str(threads), '--bowtie2',
        '--bowtie2-sensitivity-level', 'very_sensitive', genomeBam,
        estimateRspd, '--forward-prob', str(forProb), read1, index, outPrefix]
    # Process read2
    if read2:
        rsemCommand.insert(-3, '--paired-end')
        rsemCommand.insert(-2, read2)
    # Process bowtie2 arguments
    if bowtie2Path:
        rsemCommand.insert(4, bowtie2Path)
        rsemCommand.insert(4, '--bowtie2-path')
    # Concatenate and return command
    rsemCommand =  filter(None, rsemCommand)
    rsemCommand = ' '.join(rsemCommand)
    return(rsemCommand)

def tophat2Align(
        genomeIndex, transcriptIndex, outDir, read1, read2 = None,
        path = 'tophat', forProb = 0.5, mateDist = 50, mateSD = 20,
        threads = 1, sampleName = None, libraryID = None, readGroup = None,
        discordant = True, mixed = True, check = True
    ):
    ''' Function to generate and return a command to run a tophat2
    alignment of paired-end RNA seq data. Function designed with
    version 2.0.14 of Tophat2. Function takes 16 arguments:
    
    1)  genomeIndex - Full path and prefix of bowtie2 genome index.
    2)  transcriptIndex - Full path and prefix of transcript index.
    3)  outDir - Name of directory in which to store output files.
    4)  read1 - Read1 FASTQ file(s); can be file or list of files.
    5)  read2 - Read2 FASTQ file(s); can be file or list of files.
    6)  path - Full path to tophat2 executable.
    7)  forwardProb - probabiility of first read coming from the
        forward strand.
    8)  mateDist - Mean inner distance between paired reads
    9)  mateSD - Standard deviation of mateDist.
    10) threads - number of threads to use.
    11) sampleName - Sample Name
    12) librayID - ID of sequencing library
    13) readGroup - Read group.
    14) discordant - Boolean, whether to return discordant alignments.
    15) mixed - Boolean, whether to return mixed alignments.
    16) check - Boolean, whether to check for suffixes to indices.
    
    '''
    # Check for suffixes to transcriptome and genome indices
    if not isinstance(check, bool):
        raise TypeError('check argument must be bool')
    if check:
        # Check for Bowtie2 Genome indices
        suffix = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
            '.rev.2.bt2', '.fa']
        for s in suffix:
            if not os.path.isfile(genomeIndex + s):
                raise IOError('Genome index file %s not found' %(
                    genomeIndex + s
                ))
            if not os.path.isfile(transcriptIndex + s):
                raise IOError('Transcript index file %s not found' %(
                    transcriptIndex + s
                ))
    # Process forProb argument
    if not isinstance(forProb, (float, int)):
        raise TypeError('forProb argument must be numeric')
    if not 1 >= forProb >= 0:
        raise ValueError('forProb argument must be >= 0 and <= 1')
    if forProb == 1:
        forProb = 'fr-secondstrand'
    elif forProb == 0:
        forProb = 'fr-firststrand'
    else:
        forProb = 'fr-unstranded'
    # Process discordant
    if not isinstance(discordant, bool):
        raise TypeError('discordant argument must be bool')
    if discordant:
        discordant = ''
    else:
        discordant = '--no-discordant'
    # Process mixed
    if not isinstance(mixed, bool):
        raise TypeError('mixed argument must be bool')
    if mixed:
        mixed = ''
    else:
        mixed = '--no-mixed'
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if not isinstance(read1, str):
        raise TypeError('read1 argument must be string or list of strings')
    if isinstance(read2, list):
        read2 = ','.join(read2)
    if not read2 is None and not isinstance(read2, str):
        raise TypeError('read2 argument must be string or list of strings')
    # Check threads argument
    if not isinstance(threads, int):
        raise TypeError('threads argument must be integer')
    if not 32 >= threads >= 1:
        raise ValueError('threads argument must be >= 1 and <= 32')
    # Build tophat2 command
    tophatCommand = [path, '-o', outDir, '-p', str(threads),
        '--transcriptome-index', transcriptIndex, discordant,
        mixed, '--library-type', forProb, '--keep-fasta-order',
        '--b2-very-sensitive', genomeIndex, read1, read2]
    # Process mate arguments
    if not isinstance(mateDist, int):
        raise TypeError('mateDist argument must be integer')
    if not -100 <= mateDist <= 5000:
        raise ValueError('mateDist argument must be >= -100 and <= 5000')
    if not isinstance(mateSD, int):
        raise TypeError('mateSD argument must be integer')
    if not 10 <= mateSD <= 5000:
        raise ValueError('mateSD argument must be >= 10 and <= 5000')
    if read2:
        mateArgument = ['--mate-inner-dist', str(mateDist),
            '--mate-std-dev', str(mateSD)]
        tophatCommand = tophatCommand[:7] + mateArgument + tophatCommand[7:]
    # Add read group information
    if not sampleName is None:
        if not isinstance(sampleName, str):
            raise TypeError('sampleName argument must be string')
        tophatCommand.insert(-3, '--rg-sample')
        tophatCommand.insert(-3, sampleName)
    if not libraryID is None:
        if not isinstance(libraryID, str):
            raise TypeError('libraryID argument must be string')
        tophatCommand.insert(-3, '--rg-library'), 
        tophatCommand.insert(-3, libraryID)
    if not readGroup is None:
        if not isinstance(readGroup, str):
            raise TypeError('readGroup argument must be string')
        tophatCommand.insert(-3, '--rg-id',)
        tophatCommand.insert(-3, readGroup)
    # Concatenate and return command
    tophatCommand = filter(None, tophatCommand)
    tophatCommand = ' '.join(tophatCommand)
    return(tophatCommand)

