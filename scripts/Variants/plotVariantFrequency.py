'''plotChromVariantFrequency.py
    
Usage:
    plotChromVariantFrequency.py <mincov> <outfile> <snpfiles>...

'''
import collections
import os
import re
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from general_python import docopt
import numpy as np
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
if not args['<outfile>'].endswith('.png'):
    raise TypeError('Outfile must have .png suffix')

def parseInputFiles(infiles, minCov=1):
    ''' Function to extract frequency for SNPs in an input
    text file.
    
    Args:
        infile (str)- Path to input text file. File must have header with
            columns named 'Cov' and 'Freq' for the SNP coverage and SNP
            frequency, respectivly.
        
        minCov (int)- Minimum coverage of SNP required.
    
    Returns:
        snpDict - A collections default dictionary of lists. Each list
            contains tuples of the position and frequency of each variant.
            The lists are ordered by position.
    
    '''
    # Check arguments
    if not isinstance(minCov, int):
        raise TypeError('minCov must be an integer')
    if not minCov > 0:
        raise ValueError('minCov must be >=1')
    # Create ordered dictionary to store data
    freqDict = collections.OrderedDict()
    for infile in infiles:
        freqList = []
        with open(infile) as filein:
            # Extract index of key columns from header
            header = filein.next().strip().split('\t')
            covIndex = header.index('Cov')
            freqIndex = header.index('Freq')
            # Extract SNP location and frequency
            for line in filein:
                splitline = line.strip().split('\t')
                cov = splitline[covIndex]
                if cov < minCov:
                    continue
                freq = splitline[freqIndex]
                except IndexError:
                    continue
                freqList.append(float(freq))
        # Sort and return data
        freqList.sort()
        freqDict[infile] = freqList
    # Sort and return data
    return(freqDict)


def createFrequencyPlots(freqDict, outFile, bins=20, range=(0,1)):
    ''' Function to create plot of SNP frequencies across chromosomes.
    
    Args:
        snpDict - A collections.OrderedDict where the key is the file name
            name and the value is an ordered list of SNP frequencies.
        outFile - The name of the output file in which to save the data.
    
    '''
    # Extract maximum counts:
    maxCount = 0
    for freq in freqDict.values():
        counts, edges = np.histogram(freq, bins=bins, range=range)
        maxCount = max(max(counts), maxCount)
    # Loop through frequencies and create plot
    plt.figure(1, figsize=[6, 4 * len(freqDict)])
    for count, sample in enumerate(freqDict):
        # Create subplot
        plt.subplot(len(freqDict), 1, count + 1)
        plt.hist(freqDict[sample], bins=edges)
        plt.title(sample)
        # Format y-axis
        plt.ylabel('Count')
        plt.ylim([0, maxCount])
        # Format x axis
        plt.xlim([0, 1])
        plt.xlabel('Frequency')
        plt.xticks([0, 0.25, 0.5, 0.75, 1])
        plt.tight_layout()
    # Save file
    plt.savefig(outFile)

if __name__ == "__main__":
    sampleFreq = parseInputFiles(args['<snpfiles>'], int(args['<mincov>']))
    for sample in sampleFreq:
        print('{}\t{}'.format(sample, len(sampleFreq[sample])))
    createFrequencyPlots(sampleFreq, args['<outfile>'])
