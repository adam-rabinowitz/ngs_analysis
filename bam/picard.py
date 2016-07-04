def index(
        inBam, picardPath, javaPath = 'java'
    ):
    ''' Function to index BAM files using the Picard toolkit. Function
    takes three arguments:
    
    1)  inBam - Name of input BAM
    2)  picardPath - Path to Picard.jar file.
    3)  javaPath - Path to java executable.
    
    '''
    # Build and return command
    indexCommand = '%s -jar %s BuildBamIndex INPUT=%s OUTPUT=%s' %(
        javaPath,
        picardPath,
        inBam,
        inBam + '.bai'
    )
    return(indexCommand)

def markDuplicates(
        inBam, outBam, logFile, picardPath, javaPath = 'java',
        removeDuplicates = False, delete = True
    ):
    ''' Function to mark duplicates using the picard toolkit.The BAM
    output file is also indexed. Function built for picard-tools version
    1.140 and takes 7 arguments:
    
    1)  inBam - Input BAM file.
    2)  outBam - Output BAM file.
    3)  logFile - Name of file in which to same output metrics.
    4)  picardPath - Path to picard jar file.
    5)  javaPath - Path to java executable.
    6)  removeDuplicates - Boolean; whether duplicates should be
        removed from the output BAM file.
    7)  delete - Boolean whether to delete input BAM file.
    
    '''
    # Process removeDuplicates option
    if removeDuplicates:
        removeDuplicates = 'REMOVE_DUPLICATES=true'
    else:
        removeDuplicates = 'REMOVE_DUPLICATES=false'
    # Create command
    duplicateCommand = [javaPath, '-jar', picardPath, 'MarkDuplicates',
        'I=' + inBam, 'O=' + outBam, 'M=' + logFile, 'ASSUME_SORTED=true',
        'CREATE_INDEX=true', removeDuplicates]
    # merge command
    duplicateCommand = ' '.join(duplicateCommand)
    # delete input if requested
    if delete:
        duplicateCommand += ' && rm %s.ba?*' %(inBam[:-4])
    # Return command
    return(duplicateCommand)

def addReplaceReadGroup(
        inBam, outBam, name, library, platform, barcode, picardPath,
        javaPath = 'java', group = '1', delete = True
    ):
    # Create command
    rgCommand = [javaPath, '-jar', picardPath, 'AddOrReplaceReadGroups',
        'I=' + inBam, 'O=' + outBam, 'RGID=' + group, 'RGLB=' + library,
        'RGPL=' + 'platform', 'RGPU=' + barcode, 'RGSM=' + name]
    # Merge command and add deletion if required
    rgCommand = ' '.join(rgCommand)
    if delete:
        rgCommand += ' && rm {}'.format(inBam)
    # Return command
    return(rgCommand)

def sortSam(
        inBam, outBam, picardPath, javaPath = 'java', nameSort = False,
        delete = True
    ):
    # Process sort type
    if nameSort:
        sort = 'queryname'
    else:
        sort = 'coordinate'
    # Create sort command
    sortCommand = [javaPath, '-jar', picardPath, 'SortSam', 'I=' + inBam,
        'O=' + outBam, 'SORT_ORDER=' + sort]
    # Merge command and add deletion if required
    sortCommand = ' '.join(sortCommand)
    if delete:
        sortCommand += ' && rm {}'.format(inBam)
    # Return command
    return(sortCommand)

def sort_addrg_mdup(
        inBam, outBam, group, name, library, platform, barcode,
        picardPath, javaPath = 'java', removeDuplicates = False
    ):
    # Create output file names
    if not outBam.endswith('.bam'):
        raise ValueError("Output file must have '.bam' suffix")
    sortBam = outBam[:-4] + '.sort.bam'
    rgBam = outBam[:-4] + '.rg.bam'
    mdupLog = outBam[:-4] + '.dedup.log'
    # Create seperate commands
    sortCommand = sortSam(inBam = inBam, outBam = sortBam,
        picardPath = picardPath, javaPath = javaPath, nameSort = False,
        delete = True)
    rgCommand = addReplaceReadGroup(inBam = sortBam, outBam = rgBam,
        name = name, library = library, platform = platform, barcode = barcode,
        picardPath = picardPath, javaPath = javaPath, group = group,
        delete = True)
    mdupCommand = markDuplicates(
        inBam = rgBam, outBam = outBam, logFile = mdupLog,
        picardPath = picardPath, javaPath = javaPath,
        removeDuplicates = removeDuplicates, delete = True
    )
    # Concatenate commands and return
    jointCommand = ' && '.join([sortCommand, rgCommand, mdupCommand])
    return(jointCommand)

def collectInsertSizeMetrics(
        inBam, outPrefix, picardPath, javaPath = 'java'
    ):
    ''' Function to calculate insert size metrics using picard toolkit.
    Function takes 4 arguments:
    
    1)  inBam - Name of input BAM
    2)  outPrefix - Prefix of output files. A '.txt' and '.pdf' file with
        this prefix will be created
    3)  picardPath - Path to Picard.jar file.
    4)  javaPath - Path to java executable.
    
    '''
    # Build and return command
    indexCommand = '%s -jar %s CollectInsertSizeMetrics I=%s O=%s H=%s' %(
        javaPath, picardPath, inBam, outPrefix + '.txt', outPrefix + '.pdf')
    return(indexCommand)
