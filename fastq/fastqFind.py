''' Functions to identify and return names of fastq files '''
import re
import os

def findFastq(prefix, dirList, pair = True, gzip = True):
    ''' A function to identify R1 and R2 FASTQ files from directories
    using a supplied filename prefix. Function returns two list
    containing the read1 and read2 files respectively. Function has
    three arguments:
    
    1)  prefix - Prefix of the FASTQ files to be found.
    2)  dirList - A list of directories to search.
    3)  pair - A boolean indicating whether to return only paired reads. Paired
        reads must be in the same directory to be identified. Default = True.
    
    The output is a list of list where the first list is the identified
    read1 files and the second list is a list of the identified read2
    file
    
    '''
    # Check pair argument
    if not isinstance(pair, bool):
        raise ValueError("'pair' argument must be boolean")
    # Create variables to store results
    read1 = []
    read2 = []
    # Create regular expression to find files
    prefix = re.escape(prefix)
    if gzip:
        read1Pattern = re.compile(
            prefix+'.*?R1(_\\d{3}){0,1}\\.fastq\\.gz$'
        )
    else:
        read1Pattern = re.compile(
            prefix+'.*?R1(_\\d{3}){0,1}\\.fastq$'
        )
    # Loop through directories to find fastq files
    for directory in dirList:
        # Extract file names
        for (dirpath, dirnames, filenames) in os.walk(directory):
            # Loop through filenames
            for f in filenames:
                # Find files matching read1 regular expression
                if re.match(read1Pattern, f):
                    # Process pairs
                    if pair:
                        # Generate name of paired read file
                        read2File, nsub = re.subn(
                            'R1(?=(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$)',
                            'R2',
                            f
                        )
                        # Raise error if read2 filename not generated
                        if nsub != 1:
                            raise IOError('Could not generate read2 filename'\
                                ' for %s' %(f))
                        # Check for existence of paired read file
                        if read2File in filenames:
                            read1Path = os.path.join(dirpath, f)
                            read1.append(read1Path)
                            read2Path = os.path.join(dirpath, read2File)
                            read2.append(read2Path)
                        else:
                            raise IOError('Could not find read2 filename for '\
                                ' %s' %(f))
                    # Store singletons
                    else:
                        read1Path = os.path.join(dirpath, f)
                        read1.append(read1Path)
    # Sort and return data
    read1.sort()
    read2.sort()
    if pair:
        return(read1,read2)
    else:
        return(read1)

def findIlluminaFastq(prefix, dirList, pair = True):
    ''' A function to identify R1 and R2 FASTQ files from directories
    using a supplied filename prefix. Function returns two list
    containing the read1 and read2 files respectively. Function has
    three arguments:
    
    1)  prefix - Prefix of the FASTQ files to be found.
    2)  dirList - A list of directories to search.
    3)  pair - A boolean indicating whether to return only paired reads. Paired
        reads must be in the same directory to be identified. Default = True.

    The output is a list of list where the first list is the identified
    read1 files and the second list is a list of the identified read2
    file

    '''
    # Check pair argument
    if not isinstance(pair, bool):
        raise ValueError("'pair' argument must be boolean")
    # Create variables to store results
    read1 = []
    read2 = []
    # Create regular expression to find files
    read1Pattern = re.compile(
        prefix+'.*?_L\\d{3}_R1(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$'
    )
    read2Pattern = re.compile(
        prefix+'.*?_L\\d{3}_R2(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$'
    )
    # Loop through dirList to find files
    for directory in dirList:
        # Extract file names
        fileList = os.listdir(directory)
        fileList.sort()
        # Find paired reads
        if pair:
            # Loop through files
            for f in fileList:
                # Find files matching read1 regular expression
                if re.match(read1Pattern, f):
                    # Generate name of paired read file
                    read2File = re.sub(
                        'R1(?=(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$)',
                        'R2',
                        f
                    )
                    # Check for existence of paired read file
                    if read2File in fileList:
                        read1.append(directory + f)
                        read2.append(directory + read2File)
        # Or find potentially unpaired reads
        else:
            # Loop through files
            for f in fileList:
                if re.match(read1Pattern, f):
                    read1.append(directory + f)
                elif re.match(read2Pattern, f):
                    read2.append(directory + f)
    # Return data
    return(read1,read2)

def findIlluminaPrefix(directory):
    ''' Function to find all Illumina FASTQ prefixes within a specified
    directory.
    '''
    # Create output variable and regular expression
    prefix = set()
    readPattern = re.compile(
        '(.*?)_L\\d{3}_R1(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$')
    # Loop trough files and extract prefixes
    files = os.listdir(directory)
    for f in files:
        match = re.match(readPattern, f)
        if match:
            prefix.add(match.group(1))
    # Return results
    return(prefix)
