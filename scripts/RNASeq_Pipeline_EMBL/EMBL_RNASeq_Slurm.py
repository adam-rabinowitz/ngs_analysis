import collections
import os
import re
import subprocess
import time

class submitJobs(object):
   
    def __init__(self):
        self.commandDict = collections.OrderedDict()
    
    def add(
            self, command, processors = 1, memory = 4, stdout = None,
            stderr = None, modules = [], depend = [], name = None,
            partition = '1day'
        ):
        ''' command to add commands to job dictionary.
        
        Args:
            command (str)- Command to submit.
            processors (int)- Number of processors.
            memory (int)- Memory per processor in gigabytes.
            stdout (str)- Full path to stdout.
            stderr (str)- Full path to stderr.
            modules (list)- A List of modules to load.
            depend (list)- A list of job ids returned from object.
            name (str)- Name of the job.
            partition (str)- Name of partition to use
        
        Returns:
            commandNo (int)- Number of submitted command.

        '''
        # Check command
        if not isinstance(command, str):
            raise TypeError('command must be a string')
        # Check processor arguments
        if not isinstance(processors, int):
            raise TypeError('processors argument must be an integer')
        if processors < 0 or processors > 16:
            raise ValueError('processors argument must be >= 1 and <= 16')
        # Check memory arguments
        if not isinstance(memory, int):
            raise TypeError('memory must be an integer')
        if not memory > 0:
            raise ValueError('memory must be > 0')
        # Process stdout and stderr
        if stdout is None:
            if stderr is None:
                stodut = stderr = '/dev/null'
            else:
                stdout = stderr
        elif stderr is None:
            stderr = stdout
        # check stdout and stderr arguments
        stdoutDir = os.path.dirname(stdout)
        if not os.path.isdir(stdoutDir):
            raise IOError('stdout directory does not exist')
        stderrDir = os.path.dirname(stderr)
        if not os.path.isdir(stderrDir):
            raise IOError('stderr directory does not exist')
        # check modules
        if not isinstance(modules, list):
            raise IOError('modules must be a list')
        for m in modules:
            if not isinstance(m, str):
                raise IOError('modules must be a list of strings')
        # check dependencies
        if not isinstance(depend, list):
            raise IOError('dependencies must be a list')
        for d in depend:
            if d not in self.commandDict:
                raise IOError('dependencies must be present in object')
        # Check name argument
        if name is not None and not isinstance(name, str):
            raise ValueError('name must be None or a string')
        # Check queue argument
        if partition not in ('htc', '1day', '1week', '1month', 'gpu', 'bigmem'):
            raise ValueError('Unrecognised partition')
        # add command to dictionary and return command number
        commandNo = len(self.commandDict)
        self.commandDict[commandNo] = (command, processors, memory, stdout,
            stderr, modules, depend, name, partition)
        return(commandNo)
    
    def status(
            self, slurmID
        ):
        ''' Function to check status of slurm job. 
        
        Args:
            slurmID (str)- ID of slurm job.
        
        Returns:
            slurmStatus (str)- Status of slurm job
        
        '''
        # Create command to check job status
        slurmStatus = None
        sacctCommand = ['sacct', '-j', slurmID, '--format', 'JobID,State']
        # Try to get status 60 times
        for _ in range(60):
            # Extract slurm status
            sacctData = subprocess.check_output(sacctCommand)
            for line in sacctData.split('\n'):
                lineData = line.strip().split()
                if len(lineData) != 2:
                    continue
                if lineData[0] != slurmID:
                    continue
                slurmStatus = lineData[1]
                break
            # Or wait 0.5 seconds
            if slurmStatus is None:
                time.sleep(0.5)
            else:
                break
        # Return slurm status or raise an error
        if slurmStatus is None:
            raise IOError('Could not determine job status')
        return(slurmStatus)
    
    def submit(
            self, check_sub = False, verbose = False
        ):
        ''' Function to submit stored slurm jobs.
        
        Args:
            check_sub (bool)- Check status of jubs after submission.
            verbose (bool)- Return verbose output.
        
        Returns:
            slurmJobList - Either a list of tuples containing submission
                parametes or, if verbose is True, a string of job
                submission paremeters
        '''
        # Check argument
        if not isinstance(check_sub, bool):
            raise TypeError('check_sub not boolean')
        if not isinstance(verbose, bool):
            raise TypeError('verbose not boolean')
        # Create slurm ID list
        slurmJobList = []
        slurmRegx = re.compile('^\\s*Submitted batch job (\\d+)\\s*$')
        # extract commands and parameters
        for commandNo in self.commandDict:
            # Unpack command and parameters
            (command, processors, memory, stdout, stderr, modules, depend,
                name, partition) = self.commandDict[commandNo]
            # Escape $ and " characters in command
            escapeCommand = re.sub('(\$|")', '\\\\\\1', command)
            # Set required modules
            moduleCommand = 'module purge'
            if modules:
                moduleOut = ' '.join(modules)
                moduleCommand += ' && module load {}'.format(' '.join(modules))
            else:
                moduleOut = None
            completeCommand = '{} && {}'.format(moduleCommand, escapeCommand)
            # create slurm command and add node informaion
            slurmCommand = ['sbatch', '--wrap', '"'+ completeCommand + '"',
                '-c', str(processors), '-o', stdout, '-e', stderr,
                '--mem-per-cpu', str(memory) + 'G', '-p', partition]
            if name:
                slurmCommand.extend(['-J', name])
            # Add dependecies to slurm command
            dependList = []
            for d in depend:
                dependList.append(slurmJobList[d][0])
            if dependList:
                dependString = 'afterany:' + ':'.join(dependList)
                slurmCommand.extend(['-d', dependString])
                dependOut = ' '.join(dependList)
            else:
                dependOut = None
            # Try and submit complete command ten times
            slurmCommand = ' '.join(slurmCommand)
            slurmID = None
            for x in range(10):
                # Submit job and check job id returned
                slurmOut = subprocess.check_output(slurmCommand, shell=True)
                slurmOut = slurmOut.decode('utf-8')
                slurmMatch = slurmRegx.match(slurmOut)
                if slurmMatch is None:
                    time.sleep(2)
                    continue
                # Extract job id and break loop if no checks
                slurmID = slurmMatch.group(1)
                if not check_sub:
                    break
                # Check status of loop
                slurmStatus = self.status(slurmID)
                if 'FAIL' in slurmStatus.upper():
                    slurmID = None
                    time.sleep(2)
                    continue
                else:
                    break
            # Raise Error if no job submitted
            if slurmID is None:
                raise IOError('Could not submit slurm job')
            # store succesfully commoted commands
            slurmJobList.append((slurmID, name, command, processors, memory,
                stdout, stderr, moduleOut, dependOut, parition))
        # save command to output file
        return(slurmJobList)

    def output(self, slurmJobList):
        ''' Function to return output string for slurm job list.
        
        Args:
            slurmJobList (list)- A list of slurm job data generated by the
                self.submit command.
        
        Returns:
            outstring (str)- Output string to write to file.
        
        '''
        # Create and extend output string
        output = ''
        for jobTuple in slurmJobList:
            output += '{:<15}{}\n'.format('job id:', jobTuple[0])
            output += '{:<15}{}\n'.format('job name:', jobTuple[1])
            output += '{:<15}{}\n'.format('command:', jobTuple[2])
            output += '{:<15}{}\n'.format('processors:', jobTuple[3])
            output += '{:<15}{}\n'.format('memory:', jobTuple[4])
            output += '{:<15}{}\n'.format('stdout:', jobTuple[5])
            output += '{:<15}{}\n'.format('stderr:', jobTuple[6])
            output += '{:<15}{}\n'.format('modules:', jobTuple[7])
            output += '{:<15}{}\n'.format('dependencies:', jobTuple[8])
            output += '{:<15}{}\n\n'.format('dependencies', jobTuple[9])
        return(output)

def extractCommand(slurmFile):
    # Extract commands from file
    regx = re.compile(
        'job id:\s+([^\n]+)\njob name:\s+([^\n]+)\ncommand:\s+([^\n]+)\n'\
        'processors:\s+([^\n]+)\nmemory:\s+([^\n]+)\nstdout:\s+([^\n]+)\n'\
        'stderr:\s+([^\n]+)\nmodules:\s+([^\n]+)\ndependencies:\s+([^\n]+)\n'
    )
    with open(slurmFile) as infile:
        data = infile.read()
    # Check line count
    lineCount = data.count('\n')
    if lineCount == 0:
        raise IOError('slurm file is empty')
    # Extract commands and check
    commandList = regx.findall(data)
    linesPerCommand = lineCount / float(len(commandList))
    if linesPerCommand != 10:
        raise IOError('error extracting commands from slurm file')
    return(commandList)

def checkStatus(slurmID):
    # Extract state data
    sacctCommand = 'sacct -j {} --format=JobID,State'.format(slurmID)
    sacctData = subprocess.check_output(sacctCommand, shell=True)
    sacctData = sacctData.decode('utf-8')
    sacctDataList = sacctData.split()
    jobIndex = sacctDataList.index(slurmID)
    state = sacctDataList[jobIndex + 1].lower()
    # Check and return output
    if state not in ('running', 'completed', 'pending', 'failed'):
        raise ValueError('State {} not recognised'.format(state))
    return(state)

def parseSlurmFile(slurmFile):
    ''' Returns job id, job name and job status of all commands in a slurm
    output file.
    
    Args:
        slurmFile (str)- Path to slurm file.
    
    Returns:
        outList (list)- A list of three element tuples containing slurm id,
            job name and current job status.
    
    '''
    # Create output variable and extract commands
    outList = []
    commandList = extractCommand(slurmFile)
    # Populate output variable and return
    for command in commandList:
        status = checkStatus(command[0])
        outList.append((command[0], command[1], status))
    return(outList)

