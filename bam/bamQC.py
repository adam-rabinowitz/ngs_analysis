import os
from general_python import toolbox

def RNASeqC(inBam, fasta, gtf, rRNA, outDir, outPrefix, seqcPath,
    javaPath = 'java', check = True, singleEnd = False):
    ''' Function to perform RNASeqC analysis of RNA-Seq aligned BAM.
    Duplicates should be marked, but not deleted from the BAM file.
    Function takes 9 arguments.

    1)  inBam - Input BAM file with duplicates marked.
    2)  fasta - FASTA file containing genome sequence.
    3)  gtf - GTF file locating transcripts within genome.
    4)  rRNA - Interval list file locating rRNA genes.
    5)  outDir - Path to output directory.
    6)  outPrefix - Prefix to use for output files.
    7)  seqcPath - Path to RNASeQC.jar
    8)  javaPath - Path to java executable
    9)  check - Boolean indicating whther to check for required FASTA
        file indices and dictionaries.

    '''
    # Check FASTA file
    toolbox.checkArg(check, 'bool')
    if check:
        fai = fasta + '.fai'
        dict = fasta[:-2] + 'dict'
        for f in [fai, dict]:
            if not os.path.isfile(f):
                raise IOError('Genome file %s not found' %(f))
    # Create sample data string
    sampleData = "'" + outPrefix + "|" + inBam + '|' + outPrefix + "'"
    # Create command
    seqcCommand = [javaPath, '-jar', seqcPath, '-r', fasta, '-rRNA', rRNA,
        '-t', gtf, '-o', outDir, '-s', sampleData]
    # Process singleEnd Argument
    toolbox.checkArg(singleEnd, 'bool')
    if singleEnd:
        seqcCommand.append('-singleEnd')
    # Join and return command
    seqcCommand = ' '.join(seqcCommand)
    return(seqcCommand)
