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

test = collectInsertSizeMetrics(
    inBam = '/farm/scratch/rs-bio-lif/rabino01/Ascl1/ATACseq/bamFiles/SR_8.bam',
    outPrefix = '/farm/scratch/rs-bio-lif/rabino01/Ascl1/ATACseq/bamFiles/SR_8.insert',
    javaPath = '/farm/babs/redhat6/software/jdk1.7.0_40/bin/java',
    picardPath = '/farm/babs/redhat6/software/picard-tools-1.140/picard.jar'
)
print test
