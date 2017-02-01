''' Functions to identify and return names of fastq files '''
import re
import os

def findFastq(prefix, dirList, pair = True):
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
    read1Pattern = re.compile(
            prefix+'.*?R1(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$'
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

def findFastqMultiPrefix(prefixList, dirList, pair = True):
    ''' A function to identify R1 and R2 FASTQ files from directories
    using a supplied list of FASTQ prefixes.
    
    Args:
        prefixList - A tuple or list of file prefixes.
        dirList - A tuple or list of directories to search.
        pair (bool)- Whether to return only paired reads. Paired
            reads must be in the same directory. Default = True.
    
    Returns:
        fastqDict - A dictionary where the key is the prefix and the
            value lists the FASTQ files. Where pair is True value is a
            two element list of lists where the first element is read1
            files and the second element is read2 files. Where pair is
            False the value is simply a list of read1 files.
    
    '''
    # Check pair argument
    if not isinstance(prefixList, (tuple, list)):
        raise TypeError('prefixList must be tuple or list')
    if not isinstance(dirList, (tuple, list)):
        raise TypeError('dirList must be tuple or list')
    if not isinstance(pair, bool):
        raise ValueError("'pair' argument must be boolean")
    # Create dictionary to store Fastq files
    fastqDict = {}
    for prefix in prefixList:
        if pair:
            fastqDict[prefix] = [[],[]]
        else:
            fastqDict[prefix] = []
    # Create regular expressions for each prefix
    regxDict = {}
    for prefix in prefixList:
        escPrefix = re.escape(prefix)
        read1Pattern = re.compile(
            escPrefix+'.*?R1(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$')
        regxDict[prefix] = read1Pattern
    # Create regular expression for read2
    regxRead2 = re.compile('R1(?=(_\\d{3}){0,1}\\.fastq(\\.gz){0,1}$)')
    # Loop through directories to find fastq files
    for directory in dirList:
        # Extract file names
        for (dirpath, dirnames, filenames) in os.walk(directory):
                # Loop through filenames
                for f in filenames:
                    # Loop through prefixes
                    for prefix in regxDict:
                        # Find files matching read1 regular expression
                        if regxDict[prefix].match(f):
                            # Store read name
                            read1 = os.path.join(dirpath, f)
                            # Process pairs
                            if pair:
                                # Generate name of paired read file
                                read2, nsub = regxRead2.subn('R2', read1)
                                # Raise error if read2 filename not generated
                                if not nsub == 1:
                                    raise IOError('Could not generate read2 '\
                                        'for {}'.format(read1))
                                # Check for existence of paired read file
                                if not os.path.isfile(read2):
                                    raise IOError('Could not find {}'.format(
                                        read2))
                                # Store files
                                fastqDict[prefix][0].append(read1)
                                fastqDict[prefix][1].append(read2)
                            # Store singletons
                            else:
                                fastqDict[prefix].append(read1)
    # Sort and return data
    for prefix in fastqDict:
        if pair:
            for count, fastqFiles in enumerate(fastqDict[prefix]):
                fastqFiles = list(set(fastqFiles))
                fastqFiles.sort()
                fastqDict[prefix][count] = fastqFiles
        else:
            fastqFiles = fastqDict[prefix]
            fastqFile = list(set(fastqFiles))
            fastqFiles.sort()
    return(fastqDict)

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
