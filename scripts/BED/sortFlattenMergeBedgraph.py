'''sortFlattenMergeBedgraph.py
    
Usage:
    sortFlattenMergeBedgraph.py <bedgraph> <chrsizes> <bedtools> <bed>... 

'''
# Load required modules
import os
from general_python import docopt
from general_python import slurm
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
if not args['<bedgraph>'].endswith('.bedgraph'):
    raise ValueError('bedgrpah file must have .bedgraph extension')
args['<bedgraph>'] = os.path.abspath(args['<bedgraph>'])
args['<chrsizes>'] = os.path.abspath(args['<chrsizes>'])
args['<bed>'] = [os.path.abspath(x) for x in args['<bed>']]
slurmJobs = slurm.submitJobs()
# Loop through bed files and create commands to sort and flatten files
outprefix = args['<bedgraph>'][:-9]
flatFileList = []
sortFlatJobID = []
logFileList = []
for count, bed in enumerate(args['<bed>']):
    # Create and store file names
    sortFile = '{}.sort.{}.bed'.format(outprefix, count + 1)
    flatFile = '{}.flat.{}.bed'.format(outprefix, count + 1)
    flatFileList.append(flatFile)
    logFile = '{}.log.{}.txt'.format(outprefix, count + 1)
    logFileList.append(logFile)
    # Create commands
    sortCommand = 'bedtools sort -i {} > {}'.format(bed, sortFile)
    flatCommand = 'bedtools merge -i {} > {}'.format(sortFile, flatFile)
    rmCommand = 'rm {}'.format(sortFile)
    jointCommand = '{} && {} && {}'.format(sortCommand, flatCommand, rmCommand)
    # Submit and store commands
    jointJobID = slurmJobs.add(
        jointCommand, modules=[args['<bedtools>']], stdout=logFile,
        stderr=logFile
    )
    sortFlatJobID.append(jointJobID)
# Create filenames
joinFile = '{}.all.join.bed'.format(outprefix)
sortFile = '{}.all.sort.bed'.format(outprefix)
logFile = '{}.all.log.txt'.format(outprefix)
# Create commands and store
catCommand = 'cat {} > {}'.format(' '.join(flatFileList), joinFile)
rmCommand = 'rm {} {}'.format(' '.join(flatFileList), ' '.join(logFileList))
sortCommand = 'bedtools sort -i {} > {}'.format(joinFile, sortFile)
bgCommand = 'bedtools genomecov -bg -i {} -g {} > {}'.format(
    sortFile, args['<chrsizes>'], args['<bedgraph>'])
# Join and store commands
bgJointCommand = '{} && {} && {} && {}'.format(
    catCommand, rmCommand, sortCommand, bgCommand)
bgJointJobID = slurmJobs.add(
    bgJointCommand, depend=sortFlatJobID, modules=[args['<bedtools>']],
    stdout=logFile, stderr=logFile
)
# Clean up command
rmCommand = 'rm {} {} {}'.format(joinFile, sortFile, logFile)
slurmJobs.add(rmCommand, depend=[bgJointJobID])
# Submit jobs
slurmJobs.submit(verbose=True, check_sub=False)
