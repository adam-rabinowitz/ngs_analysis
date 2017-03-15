'''plotChromVariantFrequency.py
    
Usage:
    plotChromVariantFrequency.py <snpfile> <mincov> <bamfile> <outfile>

'''
import collections
import os
import re
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from general_python import docopt
# Extract and process arguments
args = docopt.docopt(__doc__,version = 'v1')
args['<snpfile>'] = os.path.abspath(args['<snpfile>'])
args['<bamfile>'] = os.path.abspath(args['<bamfile>'])
args['<outfile>'] = os.path.abspath(args['<outfile>'])
if not args['<outfile>'].endswith('.png'):
    raise TypeError('Outfile must have .png suffix')

def parseInputFile(infile, minCov=1):
    ''' Function to extract location and frequency for SNPs in an input
    tesxt file.
    
    Args:
        infile (str)- Path to input text file. File must have header with
            columns named 'SNP', 'Cov' and 'Freq' for the SNP name, SNP
            coverage and SNP frequency, respectivly.
            respectively.
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
    # Create output object and open input file
    snpDict = collections.defaultdict(list)
    with open(infile) as filein:
        # Extract index of key columns from header
        header = filein.next().strip().split('\t')
        snpIndex = header.index('SNP')
        covIndex = header.index('Cov')
        freqIndex = header.index('Freq')
        # Extract SNP location and frequency
        for line in filein:
            splitline = line.strip().split('\t')
            cov = splitline[covIndex]
            if cov < minCov:
                continue
            snp = splitline[snpIndex]
            freq = splitline[freqIndex]
            chrom, position = snp.split(':')[:2]
            snpDict[chrom].append((int(position), float(freq)))
    # Sort and return data
    for chrom in snpDict:
        snpDict[chrom] = sorted(snpDict[chrom], key=lambda x: x[0])
    return(snpDict)

def extractChromosomeLengths(bamFile):
    ''' Function to extract chromosome lengths for standard chromosomes.
    
    Args:
        bamFile - Indexed BAM file.
    
    '''
    # Create regx for matching and output list
    regx = re.compile('^(chr){0,1}(\\d+|x|y)$', re.I)
    chromList = []
    inBam = pysam.AlignmentFile(bamFile)
    # Loop through BAM files and extract lengths of matches
    for count, chrom in enumerate(inBam.references):
        if regx.match(chrom):
            chromList.append((chrom, 1, inBam.lengths[count]))
    return(chromList)

def createFrequencyPlots(snpDict, regionList, outFile):
    ''' Function to create plot of SNP frequencies across chromosomes.
    
    Args:
        snpDict - A collections.defaultdict where the key is a chromsome
            name and the values is a list of tuples. Each tuple is position
            and frequency of a SNP. The lists should be ordered by position.
        regionList - A list of tuples. Each tuple contains the chromosome,
            start and end of the regions to plot.
        outFile - The name of the output file in which to save the data.
    
    '''
    # Set figure size
    plt.figure(1, figsize=[10, 2.5 * len(regionList)])
    for count, (chrom, start, end) in enumerate(regionList):
        # Extract data and alter positions
        if len(snpDict[chrom]):
            position, frequency = zip(*snpDict[chrom])
        else:
            position, frequency = [], []
        position = [float(x) / 1e6 for x in position]
        # Create subplot
        plt.subplot(len(regionList), 1, count + 1)
        plt.scatter(position, frequency, s=10)
        plt.title('{}:{}-{}'.format(chrom, start, end))
        # Format y-axis
        plt.ylim([-0.1,1.1])
        plt.yticks([0, 0.5, 1])
        plt.ylabel('Ratio')
        # Format x axis
        plt.xlim([start / 1e6, end / 1e6])
        plt.xlabel('Position Mb')
        plt.tight_layout()
    # Save file
    plt.savefig(outFile)

if __name__ == "__main__":
    regionList = extractChromosomeLengths(args['<bamfile>'])
    snpDict = parseInputFile(args['<snpfile>'], int(args['<mincov>']))
    createFrequencyPlots(snpDict, regionList, args['<outfile>'])
