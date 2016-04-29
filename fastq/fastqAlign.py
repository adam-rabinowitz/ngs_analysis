import os
import sys
from ngs_python.bam import samtools
from general_python import toolbox

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
    toolbox.checkArg(check, 'bool')
    if check:
        suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
            '.rev.2.bt2']
        for s in suffixes:
            if not os.path.isfile(index + s):
                raise IOError('Index file %s no found' %(index + s))
    # Check and process discordant
    toolbox.checkArg(discordant, 'bool')
    if discordant:
        discordant = ''
    else:
        discordant = '--no-discordant'
    # Check mixed
    toolbox.checkArg(mixed, 'bool')
    if mixed:
        mixed = ''
    else:
        mixed = '--no-mixed'
    # Check upto argument
    toolbox.checkArg(upto, 'int', mn=1)
    # Check maximum insert argument
    toolbox.checkArg(maxInsert, 'int', mn = 1)
    # Check outut file name and generate intermediate file names
    if outFile.endswith('.sam'):
        outSam = outFile
        outBam = ''
    elif outFile.endswith('.bam'):
        outBam = outFile
        outSam = outFile[:-4] + '.sam'
    else:
        raise ValueError("'ouFile' argument must end '.sam' or '.bam'")
    # Join multiple fastq files
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
    # Create initial command
    bowtie2Command = [bowtie2Path, '--phred33', '--very-sensitive', mixed, discordant, '-p',
        str(threads), '-x', index, '-S', outSam]
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
        readGroup = 1, sampleName = None, libraryID = None, platform = None,
        markSecondary = True, check = True, samtoolsPath = 'samtools',
        memory = '2', nameSort = False
    ):
    ''' Function to generate command to perform BWA mem alignment of single
    end or paired end FASTQ files. If the supplied output file name ends with
    '.bam' then a sorted BAM file will be generated else if the file names
    ends with '.sam' a sam file is returned. Function takes the following 14
    arguments:
    
    1)  index - Path and prefix of BWA index.
    2)  outFile - Full name of output file. Must end '.sam' or '.bam'.
    3)  read1 - Read1 FASTQ file.
    4)  read2 - Read2 FASTQ file.
    5)  bwaPath - Full path to BWA exectuable.
    6)  threads - Number of threads to use.
    7)  sampleName - Name of sample to be used in header.
    8)  libraryID - Name of sample to be used in header.
    9)  readGroup - Name of read group to use in header and alignments.
    10) markSecondary - Boolean, whether to mark secondary alignments.
    11) check - Boolean, whether to check for index extensions.
    12) samtoolsPath - Full path to samtools executable.
    13) memory - Memory, in G, to use in generating BAM file.
    14) nameSort - Boolean, whether to generate a name sorted BAM file.
    
    '''
    # Check index extensions if requested
    toolbox.checkArg(check, 'bool')
    if check:
        suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        for s in suffixes:
            if not os.path.isfile(index + s):
                raise IOError('Genome index file %s no found' %(index + s))
    # Process secondary command
    toolbox.checkArg(markSecondary, 'bool')
    if markSecondary:
        markSecondary = '-M'
    else:
        markSecondary = ''
    # Check outut file name and generate intermediate file names
    if outFile.endswith('.sam'):
        outSam = outFile
        outBam = ''
    elif outFile.endswith('.bam'):
        outBam = outFile
        outSam = outFile[:-4] + '.sam'
    else:
        raise ValueError("'ouFile' argument must end '.sam' or '.bam'")
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
    # Supplement BWA command with sort command
    if outBam:
        sortCommand = samtools.sort(inFile = outSam, outFile = outBam,
            name = nameSort, memory = memory, delete = True,
            path = samtoolsPath, threads = threads)
        completeCommand = bwaCommand + ' && ' + sortCommand
    else:
        completeCommand = bwaCommand
    # Return complete command
    return(completeCommand)

def rsemBowtie2Align(
        index, outPrefix, read1, read2 = None,
        rsemPath = 'rsem-calculate-expression', bowtie2Path = '', threads = 1,
        forProb = 0.5, genomeBam = True, estimateRspd = True, check = True
    ):
    ''' Function generates and returns a command to align paired FASTQ
    files to a Bowtie2 indexed RSEM refrence transcriptome. Command
    takes 10 arguments:
    
    1)  index - Full path and prefix of transcriptome index.
    2)  outPrefix - Full path and prefix of output files.
    3)  read1 - Read1 FASTQ file(s); can be file or list of files.
    4)  read2 - Read2 FASTQ file(s); can be file or list of files.
    5)  rsemPath - Full path to rsem-calculate-expression binary
    6)  bowtie2Path - Path to directory containing bowtie2 binaries
    7)  threads - Number of threads to use in alignment
    8)  forProb - Probabiliy that the first read is derived from the
        forward strand of the transcript.
    9)  genomeBam - Boolean, whether to generate a genome BAM file.
    10) esitmateRspd - Boolean, estimate read start position distribution.
    11) check - Boolean, check if correctindex suffixes are present
    
    '''
    # Check for correct suffixes
    toolbox.checkArg(check, 'bool')
    if check:
        suffix = [
            '.grp', '.ti', '.transcripts.fa', '.seq', '.idx.fa', '.1.bt2',
            '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'
        ]
        for s in suffix:
            if not os.path.isfile(index + s):
                raise IOError('Could not find index file %s' %(
                index + s
            ))
    # Process genomeBAM argument
    toolbox.checkArg(genomeBam, 'bool')
    if genomeBam:
        genomeBam = '--output-genome-bam'
    else:
        genomeBam = ''
    # Process estimateRspd argument
    toolbox.checkArg(estimateRspd, 'bool')
    if estimateRspd:
        estimateRspd = '--estimate-rspd'
    else:
        estimateRspd = ''
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
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
        threads = 1, sampleName = '', libraryID = '', readGroup = '',
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
    toolbox.checkArg(check, 'bool')
    if check:
        # Check for Bowtie2 Genome indeces
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
    toolbox.checkArg(forProb, 'num', mn = 0, mx = 1)
    if forProb == 1:
        forProb = 'fr-secondstrand'
    elif forProb == 0:
        forProb = 'fr-firststrand'
    else:
        forProb = 'fr-unstranded'
    # Process discordant
    toolbox.checkArg(discordant, 'bool')
    if discordant:
        discordant = ''
    else:
        discordant = '--no-discordant'
    # Process mixed
    toolbox.checkArg(mixed, 'bool')
    if mixed:
        mixed = ''
    else:
        mixed = '--no-mixed'
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    toolbox.checkArg(read1, 'str')
    if isinstance(read2, list):
        read2 = ','.join(read2)
    toolbox.checkArg(read2, 'str')
    # Build tophat2 command
    tophatCommand = [path, '-o', outDir, '-p', str(threads),
        '--transcriptome-index', transcriptIndex, discordant,
        mixed, '--library-type', forProb, '--keep-fasta-order',
        '--b2-very-sensitive', genomeIndex, read1, read2]
    # Process mate arguments
    if not mateDist.startswith('$'):
        toolbox.checkArg(mateDist, 'int', mn = -100, mx = 5000)
    if not mateDist.startswith('$'):
        toolbox.checkArg(mateSD, 'int', gt = 0, mx = 5000)
    if read2:
        mateArgument = ['--mate-inner-dist', str(mateDist),
            '--mate-std-dev', str(mateSD)]
        tophatCommand = tophatCommand[:7] + mateArgument + tophatCommand[7:]
    # Add read group information
    if sampleName:
        tophatCommand.insert(-3, '--rg-sample')
        tophatCommand.insert(-3, sampleName)
    if libraryID:
        tophatCommand.insert(-3, '--rg-library'), 
        tophatCommand.insert(-3, libraryID)
    if readGroup:
        tophatCommand.insert(-3, '--rg-id',)
        tophatCommand.insert(-3, readGroup)
    # Concatenate and return command
    tophatCommand = filter(None, tophatCommand)
    tophatCommand = ' '.join(tophatCommand)
    return(tophatCommand)

