import os
import sys

def bowtie2Align(index, outSam, read1, read2 = '', path = 'bowtie2',
    threads = 1, discordant = False, mixed = False, upto = None, check = True):
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
    # check for index etensions
    if isinstance(check, bool):
        if check:
            suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                '.rev.2.bt2']
            for s in suffixes:
                if not os.path.isfile(index + s):
                    raise IOError('Index file %s no found' %(index + s))
    else:
        raise ValueError("'check' argument must be boolean") 
    # Check and process discordant
    if not isinstance(discordant,bool):
        raise IOError("'discordant' argument must be boolean")
    else:
        if discordant:
            discordant = ''
        else:
            discordant = '--no-discordant'
    # Check mixed
    if not isinstance(mixed,bool):
        raise IOError("'mixed' argument must be boolean")
    else:
        if mixed:
            mixed = ''
        else:
            mixed = '--no-mixed'
    # Join multiple fastq files
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
    # Create initial command
    bowtie2Command = [path, '--phred33', '--very-sensitive', mixed, discordant, '-p',
        str(threads), '-x', index, '-S', outSam]
    # Extend command depending on if read2 is applied
    if read2:
        bowtie2Command.extend(['-1', read1, '-2', read2])
    else:
        bowtie2Command.extend(['-U', read1])
    # Supplement additional commands
    if upto:
        bowtie2Command.extend(['-u', str(upto)])
    # Concatenate and return command
    bowtie2Command = filter(None, bowtie2Command)
    bowtie2Command = ' '.join(bowtie2Command)
    return(bowtie2Command)

def bwaMemAlign(index, outSam, read1, read2 = '', path = 'bwa',
    threads = 1, sampleName = None, libraryID = None, readGroup = None,
    markSecondary = True, check = True):
    ''' Function to generate command to perform BWA mem alignment of single
    end or paired end FASTQ files. Function takes the following 10 arguments:

    1)  index - Path and prefix of BWA index.
    2)  outSam - Full name of output SAM file.
    3)  read1 - Read1 FASTQ file.
    4)  read2 - Read2 FASTQ file.
    5)  path - Full path to BWA exectuable.
    6)  threads - Number of threads to use in alignment
    7)  sampleName - Name of sample to be used in header.
    8)  libraryID - Name of sample to be used in header.
    9)  readGroup - Name of read group to use in header and alignments.
    10) markSecondary - Boolean, whether to mark secondary alignments.
    11) check - Boolean, whether to check for correct index extensions.
    
    '''
    # Check index extensions if requested
    if isinstance(check, bool):
        if check:
            suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
            for s in suffixes:
                if not os.path.isfile(index + s):
                    raise IOError('Genome index file %s no found' %(index + s))
    else:
        raise ValueError("'check' argument must be boolean") 
    # Process secondary command
    if isinstance(markSecondary, bool):
        if markSecondary:
            markSecondary = '-M'
        else:
            markSecondary = ''
    else:
        raise ValueError("'markSecondary' argument must be boolean")
    # Create command
    bwaCommand = [path, 'mem', markSecondary ,'-t', str(threads),
        index, read1, read2]
    # Remove missing elements from coomand
    bwaCommand = filter(None,bwaCommand)
    # Add read group data
    if readGroup:
        # Create read group string
        rgString = "'@RG"
        if sampleName:
            rgString = '\\tSM:' + sampleName
        if libraryID:
            rgString += '\\tLB:' + libraryID
        rgString += '\\tID:' + readGroup + "'"
        # Add string to command
        bwaCommand.insert(2,rgString)
        bwaCommand.insert(2,'-R')
    print bwaCommand
    # Complete and return bwa command
    completeCommand = '%s > %s' %(
        ' '.join(bwaCommand),
        outSam
    )
    return(completeCommand)

def rsemBowtie2Align(index, outPrefix, read1, read2 = '',
    rsemPath = 'rsem-calculate-expression', bowtie2Path = '', threads = 1,
    forProb = 0.5, genomeBam = True, estimateRspd = True, check = True):
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
    if isinstance(check, bool):
        if check:
            suffix = ['.grp', '.ti', '.transcripts.fa', '.seq', '.chrlist',
            '.idx.fa', '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
            '.rev.2.bt2']
            for s in suffix:
                if not os.path.isfile(index + s):
                    raise IOError('Could not find index file %s' %(
                        index + s
                    ))
    else:
        raise IOError("'Check' argument must be boolean")
    # Process genomeBAM argument
    if isinstance(genomeBam, bool):
        if genomeBam:
            genomeBam = '--output-genome-bam'
        else:
            genomeBam = ''
    else:
        raise IOError("'genomeBam' argument must be boolean")
    # Process estimateRspd argument
    if isinstance(estimateRspd, bool):
        if estimateRspd:
            estimateRspd = '--estimate-rspd'
        else:
            estimateRspd = ''
    else:
        raise IOError("'estimateRspd' argument must be boolean")
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
    # Build command
    rsemCommand = [rsemPath, '-p', str(threads), '--bowtie2',
        '--bowtie2-sensitivity-level', 'very_sensitive', '--phred33',
        genomeBam, estimateRspd, '--forward-prob', str(forProb), read1, index,
        outPrefix]
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

def tophat2Align(genomeIndex, transcriptIndex, outDir, read1, read2 = '',
    path = 'tophat', forProb = 0.5, mateDist = 20, mateSD = 50, threads = 1,
    sampleName = '', libraryID = '', readGroup = '', discordant = True,
    mixed = True, check = True):
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
    if isinstance(check, bool):
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
    else:
        raise IOError("'Check' argument must be boolean")
    # Process forProb argument
    if float(forProb) == 1:
        forProb = 'fr-secondstrand'
    elif float(forProb) == 0:
        forProb = 'fr-firststrand'
    elif float(forProb) > 0 and float(forProb) < 1:
        forProb = 'fr-unstranded'
    else:
        raise IOError("'forProb' argument must be >= 0 and <= 1")
    # Process discordant
    if isinstance(discordant, bool):
        if discordant:
            discordant = ''
        else:
            discordant = '--no-discordant'
    else:
        raise IOError("'discordant' argument must be boolean")
    # Process mixed
    if isinstance(mixed, bool):
        if mixed:
            mixed = ''
        else:
            mixed = '--no-mixed'
    else:
        raise IOError("'mixed' argument must be boolean")
    # Process read arguments
    if isinstance(read1, list):
        read1 = ','.join(read1)
    if isinstance(read2, list):
        read2 = ','.join(read2)
    # Build tophat2 command
    tophat2Command = [path, '-o', outDir, '-p', str(threads),
        '--transcriptome-index', transcriptIndex, '--mate-inner-dist',
        str(mateDist), '--mate-std-dev', str(mateSD), discordant,
        mixed, '--library-type', forProb, '--keep-fasta-order',
        '--b2-very-sensitive', genomeIndex, read1, read2]
    # Add read group options to command
    if sampleName:
        tophat2Command.insert(-3, '--rg-sample')
        tophat2Command.insert(-3, sampleName)
    if libraryID:
        tophat2Command.insert(-3, '--rg-library'), 
        tophat2Command.insert(-3, libraryID)
    if readGroup:
        tophat2Command.insert(-3, '--rg-id',)
        tophat2Command.insert(-3, readGroup)
    # Concatenate and return command
    tophat2Command = filter(None, tophat2Command)
    tophat2Command = ' '.join(tophat2Command)
    return(tophat2Command)

