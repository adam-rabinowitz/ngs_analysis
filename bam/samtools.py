import subprocess
import re
import os

def sam2bam(inSam, outBam, path = 'samtools', threads = 1, 
    delete = True):
    ''' Function to convert SAM file to BAM using samtools. Function built for
    samtools version 1.2. Function takes 4 arguments:
    
    1)  inSam - Name of input SAM file.
    2)  outBam - Name of output BAM file.
    3)  samtools = Path to Samtools.
    4)  delete - Boolean indicating whether to delete input SAM file.
    
    '''
    # Generate command
    convertCommand = '%s view -bh -@ %s -o %s %s' %(
        path,
        str(threads),
        outBam,
        inSam
    )
    if delete:
        convertCommand += ' && rm %s' %(
            inSam
        )
    # Return command
    return(convertCommand)

def bam2sam(inBam, outSam, path = 'samtools', delete = True):
    ''' Function to convert BAM files to SAM files using samtools. Function
    built for samtools version 1.2. Function takes 4 arguments:
    
    1)  inSam - Name of input SAM file.
    2)  outBam - Name of output BAM file.
    3)  samtools - Path to Samtools.
    4)  delete - Boolean indicating whether to delete input SAM file.
    
    '''
    # Generate command
    convertCommand = '%s view -h -o %s %s' %(
        samtools,
        outSam,
        inBam
    )
    if delete:
        convertCommand += ' && rm %s' %(
            inSam
        )
    # Return command
    return(convertCommand)

def index(inBam, path = 'samtools'):
    ''' Function to index BAM files. Function built vor version 1.2 of
    samtools. Function takes 2 arguments:
    
    1)  inBam - Input BAM file
    2)  path - Path to samtools
    
    '''
    # Create and return command
    indexCommand = '%s index %s' %(
        path,
        inBam
    )
    return(indexCommand)

def sort(inFile, outFile, name = False, threads = 1,
    memory = '2', delete = True, path = 'samtools'):
    ''' Function to sort SAM/BAM files using samtools. Function built for
    samtoools version 1.2. If the output file is a BAM and is not sorted by
    name then the file will also be indexed. Function takes 8 arguments:

    1)  inFile - Name of input file. Must end with '.bam' or '.sam' to
        indicate file type.
    2)  outFile - Name of output file. Must end with '.bam' or '.sam' to
        indicate file type.
    3)  name - Boolean indicating whether to sort by name.
    4)  threads - Number of threads to use in sort.
    5)  memory - Memory to use in each thread.
    6)  delete- Boolean indicating whether to delete input and intermediate
        files.
    7)  path - Path to samtools executable.
    
    '''
    # Check input file
    if inFile.endswith('.bam'):
        outputCommand = ''
    elif inFile.endswith('.sam'):
        outputCommand = samtools_sam2bam(
            inSam = inFile,
            outBam = inFile[:-4] + '.bam',
            path = path,
            delete = delete,
            threads = threads
        )
        inFile = inFile[:-4] + '.bam'
    else:
        raise TypeError('File %s does not end with .bam or .sam' %(
            inFile
        ))
    # Check output file
    if outFile.endswith('.bam'):
        outFormat = 'bam'
    elif outFile.endswith('.sam'):
        outFormat = 'sam'
    else:
        raise TypeError('File %s does not end with .bam or .sam' %(
            outFile
        ))
    # Check and process name argument
    if name:
        nameSort = '-n'
    else:
        nameSort = ''
    # Generate sort command
    sortCommand = [path, 'sort', nameSort, '-m', memory, '-@', str(threads),
        '-o', outFile, '-T', outFile[:-4], '-O', outFormat , inFile]
    sortCommand = filter(None,sortCommand)
    sortCommand = ' '.join(sortCommand)
    # Add sort command to output command
    if outputCommand:
        outputCommand += ' && %s' %(
            sortCommand
        )
    else:
        outputCommand = sortCommand
    # Delete input if required
    if delete:
        outputCommand += ' && rm %s' %(
            inFile
        )
    # Index output if possible
    if not name and outFormat == 'bam':
        outputCommand += ' && %s' %(
            samtools_index(
                inBam = outFile,
                path = path
            )
        )
    # Return command
    return(outputCommand)
