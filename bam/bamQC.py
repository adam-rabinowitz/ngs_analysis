import os
from general_python import toolbox

def RNASeqC(inBam, fasta, gtf, rRNA, outDir, outPrefix, seqcPath,
    javaPath = 'java', check = True, singleEnd = False, memory=None):
    ''' Function to perform RNASeqC analysis of RNA-Seq aligned BAM.
    Duplicates should be marked, but not deleted from the BAM file.
    Function takes 9 arguments.
    
    Args:
        inBam(str) - Full path to input BAM file.
        fasta (str)- Full path to genome FASTA file.
        gtf (str)- Full path to GTF file.
        rRNA (str)- Full path to interval list file of rRNA exons.
        outDir (str)- Full path to output directory.
        outPrefix (str)- Prefix to use for output files.
        seqcPath (str)- Path to RNASeQC.jar
        javaPath (str)- Path to java executable
        check (bool)- Whether to check for required FASTA
            file indices and dictionaries.
        singleEnd (bool)- Whether library is single end.
        memory (int)- Memory to use in gigabytes.
    
    '''
    # Check FASTA file
    if not isinstance(check, bool):
        raise TypeError('check argument must be bool')
    if check:
        for f in [fasta, gtf, rRNA, fasta + '.fai', fasta[:-2] + 'dict']:
            if not os.path.isfile(f):
                raise IOError('File {} not found'.format(f))
    # Create sample data string
    sampleData = "'" + outPrefix + "|" + inBam + '|' + outPrefix + "'"
    # Create command
    seqcCommand = [javaPath, '-jar', seqcPath, '-r', fasta, '-rRNA', rRNA,
        '-t', gtf, '-o', outDir, '-s', sampleData]
    # Process singleEnd Argument
    if not isinstance(singleEnd, bool):
        raise TypeError('singleEnd argument must be bool')
    if singleEnd:
        seqcCommand.append('-singleEnd')
    # Process memory argument
    if memory:
        if not isinstance(memory, int):
            raise TypeError('memory argument must be integer')
        seqcCommand.insert(1, '-Xmx{}g'.format(memory))
    # Join and return command
    seqcCommand = ' '.join(seqcCommand)
    return(seqcCommand)
