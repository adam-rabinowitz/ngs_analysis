"""generateNormalisedMatrix.py

Usage:
    
    generateNormalisedMatrix.py <mincount> <infiles>...
        [--scope=<scope>] [--threads=<threads>]
    
    generateCountMatrix.py (-h | --help)
    
Options:
    
    --threads=<threads>  Number of threads [default: 1]
    --scope=<scope>      Scope of normalisation. Either genome ('gen'),
                         chromosome ('chr') or path to a tab delimited
                         text file outlining chromosome, start, end and
                         name of regions.
    --help               Output this message
    
"""
# Import required modules
import os
import re
import numpy as np
from ngs_python.structure import interactionMatrix, analyseInteraction
from general_python import docopt, toolbox
