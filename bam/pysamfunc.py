import pysam
import collections
import multiprocessing
from general_python import toolbox
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests

def pairDist(read1, read2):
    ''' Function calculates inner and outer distance for paired reads if
    they are determined to be concordant. Function takes two arguemnts:
    
    1)  read1 - Pysam AlignedSegment for read1.
    2)  read2 - Pysam AlignedSegment for read2.
    
    '''
    # Initialise return variable
    outer = inner = None
    # Process reads on same chromsome
    if read1.reference_id == read2.reference_id:
        # Process reads where read1 on reverse and read2 on forward strand
        if read1.is_reverse and not read2.is_reverse:
            # Check distance between potential pairs
            if (read1.reference_end >= read2.reference_end and
                read1.reference_start >= read2.reference_start):
                outer = read1.reference_end - read2.reference_start
                inner = read1.reference_start - read2.reference_end
        # Process reads where read1 on forward and read2 on reverse strand
        elif not read1.is_reverse and read2.is_reverse:
            # Check distance between potential pairs
            if (read2.reference_end >= read1.reference_end and
                read2.reference_start >= read1.reference_start):
                outer = read2.reference_end - read1.reference_start
                inner = read2.reference_start - read1.reference_end
    # Return return variable
    return(outer, inner)

def pairGenerator(bamFile, mapQ = 20):
    ''' Function extracts read pairs from  name sorted BAM files.
    Secondary alignments are ignored. Function takes two arguments:
    
    1)  bamFile - Name of input bamFile.
    2)  mapQ - Minimum mapping quality of reads.
    
    '''
    # Open bamFile and loop through
    readList = []
    currentName = ''
    inBam = pysam.AlignmentFile(bamFile, 'rb')
    while True:
        # Extract read data
        try:
            read = inBam.next()
            readName = read.query_name
        except StopIteration:
            readName = 'EndOfFile'
        # Process completed families
        if readName != currentName:
            # Count and process read-pairs
            if len(readList) == 2:
                read1, read2 = readList
                if read1.is_read1 and read2.is_read2:
                    yield(read1, read2)
                elif read1.is_read2 and read2.is_read1:
                    yield(read2, read1)
            # Reset read list and current name
            currentName = readName
            readList = []
        # Break loop at end of BAM file
        if readName == 'EndOfFile':
            break
        # Skip secondary alignments
        elif (read.flag & 256):
            continue
        # Skip reads below supplied mapping quality
        if read.mapping_quality <= mapQ:
            continue
        # Else append read to read list
        else:
            readList.append(read)

def baseCalls(
        sequence, quality, positions, cigartuple, groupDeletions = True
    ):
    # Remove clipping from cigar and extract cigar types
    cigartuple = [x for x in cigartuple if x[0] < 4]
    # Check cigar types
    cigartypes = [x[0] for x in cigartuple]
    # Process indels
    if 1 in cigartypes or 2 in cigartypes:
        # Trim clipping and indels from start of cigartuple
        while cigartuple:
            # Remove indels not flanked by mapped sequence
            if cigartuple[0][0] in [1,2]:
                sequence = sequence[cigartuple[0][1]:]
                quality = quality[cigartuple[0][1]:]
                cigartuple = cigartuple[1:]
            # Stop trimming
            else:
                break
        # Remove clipping and indels from end of read
        while cigartuple:
            # Remove indels not flanked by mapped sequence
            if cigartuple[-1][0] in [1,2]:
                sequence = sequence[:-cigartuple[-1][1]]
                quality = quality[:-cigartuple[-1][1]]
                cigartuple = cigartuple[:-1]
            # Stop trimming
            else:
                break
        # Initialise variables to extract indel data
        location = 0
        insertion = []
        deletion = []
        # Find position and length of indels
        if cigartuple:
            for element in cigartuple:
                # Process mapped elements of cigar
                if element[0] == 0:
                    location += element[1]
                # Process insertion elements of cigar
                elif element[0] == 1:
                    insertion.append((location - 1, location + element[1] ))
                # Process deletion elements of cigar
                elif element[0] == 2:
                    deletion.append((location, element[1] ))
        # Process insertions
        for i in insertion:
            # Concatenate insertions into single list element
            sequence[i[0]:i[1]] = ["".join(sequence[i[0]:i[1]])]
            # Alter quality to mean of bases
            quality[i[0]:i[1]] = [ sum(quality[i[0]:i[1]]) / 
                len(quality[i[0]:i[1]]) ]
        # Reverse and loop through deletion
        deletion.reverse()
        for d in deletion:
            # Alter sequence to include "-" to signify deletion
            sequence = sequence[:d[0]] + ['-'] * d[1] + sequence[d[0]:]
            # Alter quality to mean of flank
            dQual = (quality[d[0] - 1] + quality[d[0]]) / 2
            quality = quality[:d[0]] + [dQual] * d[1] + quality[d[0]:]
            # Alter positions
            dPos = range(positions[d[0] - 1] + 1,
                positions[d[0] - 1] + d[1] + 1)
            positions = positions[:d[0]] + dPos + positions[d[0]:]
    # Raise error if reference skipped
    elif 3 in cigartypes:
        raise IOError('Irregular cigar found: %s' %(cigartuple))
    # Create and return output
    if len(sequence) == len(quality) and len(quality) == len(positions):
        positionDict = dict(zip(positions, zip(sequence, quality)))
        return(positionDict)
    else:
        raise IOError('Could not ascribe sequence to positions')
    # Return data

def extractPosition(
       openBam, chrom, position, minMapQ = 20, minBaseQ = 20
    ):
    # check arguments
    toolbox.checkArg(chrom, 'str')
    toolbox.checkArg(position, 'int', mn = 1)
    toolbox.checkArg(minMapQ, 'int', mn = 0)
    toolbox.checkArg(minBaseQ, 'int', mn = 0)
    toolbox.checkArg(groupDeletions,  'bool')
    # Set variables for mapping
    mapQuality = []
    baseCounts = {}
    # Loop through reads covering BAM
    for read in openBam.fetch(chrom, position-1, position):
        # Store mapping quality and skip reads with low values
        mapQuality.append(read.mapping_quality)
        if mapQuality[:-1] < minMapQ:
            continue
        # Extract base calls for read
        baseDict = baseCalls(
            sequence = list(read.query_alignment_sequence),
            quality = list(read.query_alignment_qualities),
            positions = list(read.get_reference_positions()),
            cigartuple = read.cigartuples
        )
        # Extract reads for base of interest or skip
        try:
            base, quality = baseDict[position - 1]
        except KeyError:
            continue
        # Skip bases of poor quality
        if quality < minBaseQ:
            continue
        # Add base to base dictionary
        if base in baseCounts:
            baseCounts[base][0] += 1
        else:
            baseCounts[base] = [1, 0]
        # Add forward strand count to base dictionary
        if not read.is_reverse:
            baseCounts[base][1] += 1
    # Calculate and return results
    if len(mapQuality) > 0:
        meanMap = sum(mapQuality) / len(mapQuality)
    else:
        meanMap = 0
    return(baseCounts, meanMap)

def extractVariantCounts(
        variantList, bamFile, minMapQ = 20, minBaseQ = 20
    ):
    ''' Function calculates the frequency at which specified nucleotides
    are found at specific chromosomal regions. Function takes 6 arguments:
    
    1)  variantList - a tuple/list of tuples/lists that contain four
        elements; the chromosome, position, reference, and variant e.g
        [('chr1', 1, 'A', 'T'), ('chr2', 100, 'C', 'G')].
    2)  bamFile - Full path to BAM file
    3)  minMapQ - Minimum mapping quality of read to extract base.
    4)  MinBaseQ - Minimum base quality to extract base.
    
    Function returns a pandas dataframe containing the following five columns:
    
    1)  varcount - Count of the variant calls.
    2)  refcount - Count of the reference calls.
    3)  varfor - Frequency of reference reads on forward strand.
    4)  reffor - Frequency of variant reads on forward strand.
    5)  mapqual - Mean mapping score of ALL reads spanning the postion.
    
    '''
    # Create output dataframe
    variantNames = [':'.join(map(str, x)) for x in variantList]
    outData = pd.DataFrame(
        columns = ['varcount', 'refcount', 'varfor', 'reffor', 'mapqual'],
        index = variantNames
    )
    # Open bamFile
    bam = pysam.AlignmentFile(bamFile)
    # Loop through variants and extract counts
    for name, (chrom, position, reference, variant) in zip(variantNames,
        variantList):
        # Extract base data
        baseCounts, mapqual = extractPosition(
            openBam = bam,
            chrom = chrom,
            position = position,
            minMapQ = minMapQ,
            minBaseQ = minBaseQ
        )
        # Count total bases
        allCount = sum([x[0] for x in baseCounts.values()])
        # Extract reference counts
        if reference in baseCounts:
            refcount, forcount = baseCounts[reference]
            reffor = round(forcount / float(allCount), 4)
        else:
            refcount = 0
            reffor = 0
        # Extract frequency
        if variant in baseCounts:
            varcount, forcount = baseCounts[variant]
            varfor = round(forCount / float(allCount), 4)
        else:
            varCount = 0
            forFreq = 0
        # Add data to output list
        outData.loc[name] = [varcount, refcount, varfor, reffor mapqual]
    # Close bam and return results
    bam.close()
    return(outData)

def calculateVariantMetrics(
        variantList, bamList, sampleNames, minMapQ = 20, minBaseQ = 20
    ):
    ''' Function calculates the frequency at which specified nucleotides
    are found at specific chromosomal regions. Function takes 6 arguments:
    
    1)  variantList - a tuple/list of tuples/lists that contain three
        elements. These three elements are the chromosome, position and
        variant of interest e.g [('chr1', 1, 'A'), ('chr2', 100, 'C')].
    2)  bamList - Full path to BAM files. It is presumed that the first
        file in the list is the reference BAM file.
    3)  sampleNames - A list/tuple of samples names. Must be one for every
        BAM files
    4)  minMapQ - Minimum mapping quality of read to extract base.
    5)  MinBaseQ - Minimum base quality to extract base.
    
    Function returns a list of tuples where each tuple contains the following
    four elements:
    
    1)  Count of the variant at the described position.
    2)  Count of the variant on the forward strand.
    3)  Frequency of base forward strand reads.
    4)  Mean mapping score of ALL reads spanning the described postion.
        This value is calculated from reads with mapping qualities and base
        qualities below the desired threshold.
    
    '''
    # Check supplied Names
    if not isinstance(sampleNames, (list, tuple)):
        raise IOError('sampleNasmes must be a list or a tuple')
    if len(sampleNames) != len(bamList):
        raise IOError('Must be a sample name for each BAM')
    # Extract variants for each BAM
    outputData = pd.DataFrame()
    for number, (name, bam) in enumerate(zip(sampleNames, bamList)):
        # Extract sample data
        sampleData = extractVariantCounts(
            variantList = variantList,
            bamFile = bam,
            minMapQ = minMapQ,
            minBaseQ = minBaseQ
        )
        # Calculate adjusted p-values
        if number:
            # Extract data
            testVar = sampleData['varcount']
            testRef = sampleData['refcount']
            # Calculate p-value
            pvalue = []
            for x in zip(zip(refVar, refRef), zip(testVar, testRef)):
                pvalue.append(fisher_exact(x)[1])
            sampleData['pvalue'] = pvalue
            # Calculate FDR
            sampleData['padj'] = multipletests(pvalue, method = 'fdr_bh')[1]
        # Extract count for reference
        else:
            refVar = sampleData['varcount']
            refRef = sampleData['refcount']
        # Rename tables columns and append to output
        sampleData.columns = [name + '_' + x for x in sampleData.columns]
        outputData = pd.concat([outputData, sampleData], axis = 1)
    # Return data
    return(outputData)

