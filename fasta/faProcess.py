from Bio import SeqIO
from general_python import toolbox
import re

def extract(fasta, chrom, start = None, end = None):
    ''' Function to extract sequence from a FASTA file which is then
    returned as a string. Function takes four arguments:
    
    1)  fasta - Input fasta file.
    2)  chrom - Chromosome name.
    3)  start - First base of sequence to extract.
    4)  end - Last base of sequence to extract.
    
    '''
    # Check arguments
    toolbox.checkArg(chrom, 'str')
    toolbox.checkArg(start, 'int')
    toolbox.checkArg(end, 'int')
    # Create iterator and loop through sequences
    for sequence in SeqIO.parse(open(fasta), 'fasta'):
        # Find chromosome within fasta file
        if sequence.id == chrom:
            # Create start and end values if not supplied
            if start is None:
                start = 1
            if end is None:
                end = len(sequence)
            # Check end is not greater than sequence length
            if end > len(sequence):
                raise ValueError('Interval extends beyond chromosome')
            # Check end is greater than or equal to start
            if end < start:
                raise ValueError('End preceeds start')
            # Extract and return sequence
            return(sequence.seq[start-1:end].tostring())
    # Raise error if chromosome not found
    raise ValueError('Chromosome not found in FASTA file')

def replace(insert, target, start, end):
    ''' Function replaces a portion of a sequence with an alternative
    sequence. Function takes four
    arguments:
    
    1)  insert - Sequence to insert into target sequence.
    2)  target - Target sequence.
    3)  start - Start of target sequence to be replaced.
    4)  end - End of target sequence to be replaces.

    Function return the novel chimeric sequence as well as the excised
    portion of the target sequence
    
    '''
    # Adjust start to account for python zero index
    start -= 1
    # Split target
    left = target[ : start]
    replace = target[start : end]
    right = target[end : ]
    # Insert reference
    chimera = left + insert + right
    return(chimera, replace)

def wrap(sequence, width = 60):
    ''' Function wraps string into lines of a fixed width. Function
    takes two arguments:
    
    1)  sequence
    2)  width - Maximum line width excluding newline character
    
    '''
    # Remove whitespace from sequence
    sequence = re.sub('\s+','',sequence)
    # Add line breaks and return
    regx = '([^\s]{%s})' %(width)
    sequence = re.sub(regx, '\\1\n', sequence)
    return(sequence.strip())
