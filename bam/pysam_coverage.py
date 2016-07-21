import pysam
import numpy as np
import pandas as pd
from general_python import toolbox
import multiprocessing 

def create_bins(sizes, width, fixed = True, overlap = 0, index = 0):
    # Check index argument
    if not index in (0, 1):
        raise ValueError('index argument must be 0 or 1')
    # Create output object and calculate span
    bins = pd.DataFrame(columns = ['chrom','start','end'])
    span = width - overlap
    count = 0
    # Loop through chromsome and generate bin data frame
    for chrom, length in sizes:
        length = float(length)
        # Generate bin start and end for equal width bins
        if fixed:
            # Calculate end of first bin
            remainder = (length - width) % span
            end = np.floor(remainder/2) + width
            # Calculate start and end of all bins
            binEnd = np.arange(end, length + 1, span, dtype = np.uint32)
            binStart = binEnd - width
        # Generate bin start and end sites for unequal width bins
        else:
            # Calculate number of bins
            binNo = np.ceil((length - width) / span) + 1
            # Calculate excess number of elements within bins
            excess = (width + (span * (binNo - 1))) - length
            # Calculate bin size 
            largeBin = int(width - np.floor(excess / binNo))
            smallBin = int(width - np.ceil(excess / binNo))
            # Calculate bin numbers
            smallBinNo = int(excess % binNo)
            largeBinNo = int(binNo - smallBinNo)
            # Generate bin width
            binWidth = np.array([largeBin] * largeBinNo + 
                [smallBin] * smallBinNo, dtype = np.uint32)
            # Generate bins
            binEnd = np.cumsum(binWidth - overlap) + overlap
            binStart = binEnd - binWidth
        # Create chromosome dataframe
        chromBins = pd.DataFrame(
            {'chrom':chrom, 'start':binStart, 'end':binEnd},
            index = np.arange(count, count + len(binStart))
        )
        if index == 1:
            chromBins['start'] += 1
        # Adjust bin count and store dataframe
        count += len(binStart)
        bins = pd.concat((bins, chromBins), axis=0)
    # Return dataframe
    bins = bins[['chrom','start','end']]
    return(bins)

def create_bins_bam(
    bam, width, fixed = True, overlap = 0, index = 1
):
    # Create chromsome list
    chrList = []
    bam = pysam.AlignmentFile(bam)
    for chrTuple in zip(bam.references, bam.lengths):
        chrList.append(chrTuple)
    bam.close()
    # Create and return bed file
    bins = create_bins(sizes = chrList, width = width, fixed = fixed,
        overlap = overlap, index = index)
    return(bins)
    
def count_reads_region(
        bam, chrom, start, end, overlap = 'any', mapq = 0, index = 0
    ):
    # Check index and adjust start coordinates
    if not index in (0, 1):
        raise ValueError('Index must be 0 or 1')
    start -= index
    # Initialize counter
    counter = 0
    # Loop through counter
    for read in bam.fetch(chrom, start, end):
        # Skip unmapped and poorly mapped reads
        if read.is_unmapped:
            continue
        if read.mapping_quality < mapq:
            continue
        # Count if overlap is 'any'
        if overlap == 'any':
            counter += 1
        # Process alternative overlaps
        else:
            # Extract position
            if overlap == 'start':
                if read.is_reverse:
                    position = read.reference_end
                else:
                    position = read.reference_start
            elif overlap == 'end':
                if read.is_reverse:
                    position = read.reference_start
                else:
                    position = read.reference_end
            else:
                raise ValueError('Unrecognised overlap argument')
            # Count if position is in desired interval
            if start <= position <= end:
                counter += 1
    # Return counter
    return(counter)

def count_reads_bam(
        bam, bed, overlap = 'any', index = 0, mapq = 0, pipe = None
    ):
    # check bed indices are unique
    if len(bed.index) > len(set(bed.index)):
        raise ValueError('bed indices must be unique')
    # Open bam file and create index
    bam = pysam.AlignmentFile(bam)
    output = pd.Series(index = bed.index, dtype = np.int32)
    # Loop through bed and extract read counts
    for dfindex, chrom, start, end in bed.itertuples():
        output[dfindex] = count_reads_region(bam = bam, chrom = chrom,
            start = start, end = end, overlap = overlap, mapq = mapq, 
            index = index)
    # Close bam and return values
    bam.close()
    if pipe:
        pipe[0].close()
        pipe[1].send(output)
        pipe[1].close()
    else:
        return(output)

def count_reads_multibam(
        bamList, bed, overlap = 'any', index = 0, mapq = 0
    ):
    # Create output dataframe
    output = pd.DataFrame(index = bed.index, columns = bamList,
        dtype = np.int32)
    # Create and store process for each bam
    processList = []
    for bam in bamList:
        pipe = multiprocessing.Pipe()
        process = multiprocessing.Process(
            target = count_reads_bam,
            args = (bam, bed, overlap, index, mapq, pipe)
        )
        process.start()
        pipe[1].close()
        processList.append((bam, process, pipe[0]))
    # Extract data for each bam
    for bam, process, pipe in processList:
        data = pipe.recv()
        pipe.close()
        process.join()
        output[bam] = data
    # Return output
    return(output)

bamList = [
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557175_sort.bam',
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557176_sort.bam',
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1557178_sort.bam',
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1559300_sort.bam',
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1564297_sort.bam',
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/SRR1564298_sort.bam',
]

bed = create_bins_bam(
    bam = bamList[0],
    width = 1000,
    fixed = True,
    index = 1
)
data = count_reads_multibam(
    bamList = bamList,
    bed = bed,
    overlap = 'end',
    mapq = 0,
    index = 1
)
data = pd.concat((bed,data), axis = 1)
data.to_csv(
    '/farm/scratch/rs-bio-lif/rabino01/Yasu/ChIP/Sutani/binCount.txt',
    sep = '\t',
    header = True,
    index = False,
    float_format = '%.0f'
)
