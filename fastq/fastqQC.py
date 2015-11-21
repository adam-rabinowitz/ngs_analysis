import re

def fastQC(inFile, outDir, path = 'fastqc'):
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
    fastqcCommand = '%s --extract -q -o %s %s; rm %s %s' %(
        path,
        outDir,
        inFile,
        outDir + name + '_fastqc.html',
        outDir + name + '_fastqc.zip'
    )
    # Execute or return command
    return(fastqcCommand)
