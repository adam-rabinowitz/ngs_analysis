import subprocess
import re
import os
from general_python import toolbox

def sam2bam(
        inSam, outBam, path = 'samtools', threads = 1, delete = True
    ):
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

def bam2sam(
        inBam, outSam, path = 'samtools', delete = True
    ):
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

def index(
        inBam, path = 'samtools'
    ):
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

def sort(
        inFile, outFile, name = False, threads = 1, memory = 2, delete = True,
        path = 'samtools'
    ):
    ''' Function to sort SAM/BAM files using samtools. Function built for
    samtoools version 1.2. If the output file is a BAM and is not sorted by
    name then the file will also be indexed. Function takes 8 arguments:

    1)  inFile - Name of input file. Must end with '.bam' or '.sam' to
        indicate file type.
    2)  outFile - Name of output file. Must end with '.bam'.
    3)  name - Boolean indicating whether to sort by name.
    4)  threads - Number of threads to use in sort.
    5)  memory - Memory, in gigabytes, to use in each thread.
    6)  delete- Boolean indicating whether to delete input file.
    7)  path - Path to samtools executable.
    
    '''
    # Check input file
    if inFile.endswith('.bam'):
        tempBam = ''
        outputCommand = ''
        sortFile = inFile
    elif inFile.endswith('.sam'):
        tempBam = inFile[:-4] + '.temp.bam'
        outputCommand = sam2bam(
            inSam = inFile,
            outBam = tempBam,
            path = path,
            delete = False,
            threads = threads
        )
        sortFile = tempBam
    else:
        raise TypeError('Input file suffix must be .bam/.sam')
    # Check output file
    if outFile.endswith('.sam'):
        outFormat = 'sam'
    elif outFile.endswith('.bam'):
        outFormat = 'bam'
    else:
        raise TypeError("Output file %s does not end '.sam/.bam'" %(outFile))
    # Check and process name argument
    if name:
        nameSort = '-n'
    else:
        nameSort = ''
    # Process memory argument
    memory = str(memory) + 'G'
    # Generate sort command
    sortCommand = [path, 'sort', nameSort, '-m', memory, '-@', str(threads),
        '-o', outFile, '-T', outFile[:-4], '-O', outFormat , sortFile]
    sortCommand = filter(None, sortCommand)
    sortCommand = ' '.join(sortCommand)
    # Add sort command to output command
    if outputCommand:
        outputCommand += ' && %s' %(
            sortCommand
        )
    else:
        outputCommand = sortCommand
    # Delete intermediate and temporary files if required
    if delete and tempBam:
        outputCommand += ' && rm %s %s' %(inFile, tempBam)
    elif delete:
        outputCommand += ' && rm %s' %(inFile)
    elif tempBam:
        outputCommand += ' && rm %s' %(tempBam)
    # Index output if possible
    if not name and outFormat == 'bam':
        outputCommand += ' && %s' %(
            index(
                inBam = outFile,
                path = path
            )
        )
    # Return command
    return(outputCommand)
    
def mpileup(
        inBam, outFile, reference = '', minMapQ = 20, minBaseQ = 20,
        countOrphans = False, countDup = False, disableBAQ = True,
        maxDepth = 10000, minFilter = 0, path = 'samtools'
    ):
    # Check numeric arguments
    if not isinstance(minMapQ, int):
        raise TypeError('minMapQ must be integer')
    if minMapQ < 2:
        raise ValueError('minMapQ must be >=2')
    if not isinstance(minBaseQ, int):
        raise TypeError('minBaseQ must be integer')
    if minBaseQ < 0:
        raise ValueError('minMapQ must be >=0')
    if not isinstance(countOrphans, bool):
        raise TypeError('countOrphans must be boolean')
    if not isinstance(countDup, bool):
        raise TypeError('countDup must be boolean')
    if not isinstance(disableBAQ, bool):
        raise TypeError('disableBAQ must be boolean')
    if not isinstance(maxDepth, int):
        raise TypeError('maxDepth must be integer')
    if maxDepth < 1:
        raise ValueError('maxDepth must be >=1')
    if not isinstance(minFilter, int):
        raise TypeError('minFilter must be integer')
    if minFilter < 0:
        raise ValueError('minFilter must be >=0')
    # Create initial command
    command = [path, 'mpileup', '-q', str(minMapQ), '-Q', str(minBaseQ), '-d',
        str(maxDepth)]
    # Add otptions
    if countOrphans:
        command.append('-A')
    if disableBAQ:
        command.append('-B')
    if reference:
        command.extend(['-f', reference])
    if not countDup:
        command.extend(['--ff', '1024'])
    # Process filter
    if minFilter:
        command.append(inBam)
        awk = 'awk \'BEGIN{FS = "\\t"};{if ($4 >= %s) print $0}\'' %(minFilter)
        finalCommand = '%s | %s > %s' %(' '.join(command), awk, outFile)
    else:
        command.extend(['-o', outFile, inBam])
        finalCommand = '%s > %s' %(' '.join(command), outFile)
    # Retutn command
    return(finalCommand)

def merge(
        inBamList, outBam, delete = False, path = 'samtools'
    ):
    # Check arguments
    if not isinstance(inBamList, list):
        raise TypeError('inBamList must be list')
    if not isinstance(delete, bool):
        raise TypeError('delete must be boolean')
    # Create and return command
    command = '{} merge -f {} {}'.format(path, outBam, ' '.join(inBamList))
    if delete:
        command += ' && rm {}'.format()
    return(command)
