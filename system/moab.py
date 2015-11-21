# Import required modules
import subprocess
import time
import re

def submitjob(command, processor = 1, stdout = '/dev/null',
    stderr = '/dev/null', dependency = None):
    ''' This function submits jobs to the farm using the MOAB software
    package and return the MOAB job identifier. The function takes 5
    arguments:
    
    1)  command - The system command to submit to the farm.
    2)  processor - Number of processors to run the job on (Default = 1).
    3)  stdout - File in which to save STDOUT stream (Default =
        /dev/null).
    4)  stderr - File in which to save STERR stream (Default =
        /dev/null).
    5)  dependency - MOAB job IDs of jobs that have to be sucessfully
        completed before processing current command (Default = None).
    
    IF 'stdout' and 'stderr' are identical then the two streams are
    combined and saved in a single file. Function returns the MOAB ID
    of the submitted job.
    '''
    # Modify input command
    if isinstance(command, list):
        command = ' '.join(command)
    # Create msub command
    msubCommand = ['/opt/moab/bin/msub', '-l']
    # Add node information
    nodeCommand = 'nodes=1:babs:ppn=%s' %(processor)
    msubCommand.append(nodeCommand)
    # Add output information
    if stdout == stderr:
        msubCommand.extend(['-j', 'oe', '-o', stdout])
    else:
        msubCommand.extend(['-o', stdout, '-e', stderr])
    # Add dependecy
    if dependency:
        if isinstance(dependency, list):
            dependency = 'x=depend:afterok:%s' %(':'.join(dependency))
        elif isinstance(dependency, str):
            dependency = 'x=depend:afterok:%s' %(dependency)
        msubCommand.extend(['-W', dependency])
    # Create output variable for function
    moabID = None
    # Try submit the job ten times
    for _ in range(10):
        # Create msub process
        msubProcess = subprocess.Popen(
            msubCommand,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        # Submit command
        moab = msubProcess.communicate(input=command)[0]
        # Search for returned Moab ID
        moabMatch = re.match('^\s+(Moab\\.\d+)\s+$',moab)
        # Check that Moab ID has been returned or wait and repeat
        if moabMatch:
            moabID = moabMatch.group(1)
            break
        else:
            time.sleep(10)
    # Stop process and return moabID
    return(moabID)
