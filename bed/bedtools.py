from general_python import toolbox

def sortBed(inBed, outBed, path = 'bedtools', delete = True):
    ''' Function to sort bed files, first by cheomsome and then by start
    position. Function takes four arguments
    
    1)  inBed - Input bed file.
    2)  outBed - Output bed file.
    3)  path - path to bedtools exectuable.
    4)  delete - Whether to delete input bed file
    
    '''
    # Check arguments
    toolbox.checkArg(delete, 'bool')
    # Create sort command:
    sortCommand = '%s sort -i %s > %s' %(path, inBed, outBed)
    # Append deletion command
    if delete == True:
        sortCommand += ' && rm %s' %(inBed)
    # Return command
    return(sortCommand)

def bed2bedGraph(inBed, outBG, chrFile, path = 'bedtools', delete = False):
    ''' Function to create bedgraph from sorted bed file. Function takes
    five arguments:
    
    1)  inBed - Input bed file.
    2)  outBG - Output bedgraph file.
    3)  chrFile - A tab delimited text file of chromosome name and sizes.
    3)  path - Path to bedtools executable.
    4)  delete - Boolean, whether to delete input BED file
    
    '''
    # Check arguments
    toolbox.checkArg(delete, 'bool')
    # Create bedgraph command
    bgCommand = '%s genomecov -bg -i %s -g %s > %s' %(path, inBed, chrFile,
        outBG)
    # Append deletion command
    if delete == True:
        bgCommand += ' && rm %s' %(inBed)
    # Return command
    return(bgCommand)
