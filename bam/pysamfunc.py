import pysam
import collections
import multiprocessing
from general_python import toolbox
from ngs_python.variant import annovar
import pandas as pd
from scipy.stats import fisher_exact
import re

def extract_fasta(fasta, chrom, start, end):
    ''' Function to extract sequence from FASTA file using the pysam
    module.
    
    Args:
        fasta: string of full path to faidx indexed FASTA file or open
            pysam.FastaFile.
        chrom (str): name of chromosome.
        start (int): start of sequence to extract (1-based index).
        end (int): end of sequence to extract (1-based index).
        
    Returns:
        str: a string of the reference sequence.
    
    Raises:
        ValueError: If desired interval not contained on chromsome.
    
    '''
    # Check arguments
    toolbox.check_var(chrom, 'str')
    toolbox.check_var(start, 'int', mn = 1)
    toolbox.check_var(end, 'int', mn = start)
    # Open FASTA if string supplied
    if isinstance(fasta, str):
        fasta = pysam.FastaFile(fasta)
    # Extract chromosome length and check end value
    chromLength = fasta.get_reference_length(chrom)
    if end > chromLength:
        raise ValueError('Interval extends beyond chromosome')
    # Extract and return sequence
    seq = fasta.fetch(chrom, start - 1, end)
    return(seq)
    
def indel_homopolymer(sequence, position, ref, var, insertion):
    # Check that insertion argument is boolean
    toolbox.check_var(sequence, 'str')
    toolbox.check_var(position, 'int')
    toolbox.check_var(ref, 'str')
    toolbox.check_var(var, 'str')
    toolbox.check_var(insertion, 'bool')
    # Covert sequences to uppercase
    sequence = sequence.upper()
    ref = ref.upper()
    var = var.upper()
    # Check that first base of indel is centre of reference sequence
    if sequence[position] != ref:
        raise ValueError('Position must be first base of indel (1-based index)')
    # Check that deletions are contained within the reference
    if not insertion and len(sequence) - position - len(var) < -1:
        raise ValueError('Deletion must be contained within reference')
    # Find homopolymer for insertion
    if insertion:
        # Extract insertion sequence and find any sub monomers
        indelSeq = var[1:]
        indelMonomer, monomerNo = toolbox.find_monomer(indelSeq)
        # Add insertion to reference
        inserted = sequence[:position] + indelMonomer + sequence[position:]
        # Create regex to find homopolymer
        regex = '^.{0,%s}?((%s)+).{0,%s}?$' %(
            position,
            indelMonomer,
            len(sequence) - position
        )
        # Find homopolymer and extract length
        homo = re.match(regex, inserted).group(1)
        homoLengthRef = (len(homo) / len(indelMonomer)) - 1
        homoLengthVar = homoLengthRef + monomerNo
    # Find homopoly
    else:
        # Extract deleted sequence
        indelSeq = sequence[position:][:len(var)]
        indelMonomer, monomerNo = toolbox.find_monomer(indelSeq)
        # Create regex to find homopolymer
        regex = '^.{0,%s}?((%s)+).{0,%s}?$' %(
            position,
            indelMonomer,
            len(sequence) - position - len(indelSeq)
        )
        # Find homopolymer and extract length
        homo = re.match(regex, sequence).group(1)
        homoLengthRef = len(homo) / len(indelMonomer)
        homoLengthVar = homoLengthRef - monomerNo
    # Return data
    return(indelMonomer, homoLengthRef, homoLengthVar)

def comp_annotate(fasta, varList):
    # Create output dataframe
    varNames = [':'.join(map(str,x)) for x in varList]
    outputDF = pd.DataFrame(index = varNames, columns = ['low_comp'])
    # Open FASTA file and create chromosome dictionary
    fastaFile = pysam.FastaFile(fasta)
    chromSize = {}
    for chrom in fastaFile.references:
        chromSize[chrom] = fastaFile.get_reference_length(chrom)
    # Loop through indel list and extract annotation
    for name, (chrom, pos, ref, var) in zip(varNames, varList):
        # Check variant is on chromosome
        if not 1 <= pos <= chromSize[chrom]:
            raise ValueError('position %s chromosome %s' %(pos, chrom))
        # Generate start and end intervals
        pos = int(pos)
        if not 1 <= pos <= chromSize[chrom]:
            raise ValueError('position {} chromosome {}'.format(pos, chrom))
        if 1 < pos < chromSize[chrom]:
            start = pos - 1
            end = pos + 1
            varIndex = 1
        elif pos == 1:
            start = pos
            end = pos + 1
            varIndex = 0
        elif pos == chromSize[chrom]:
            start = pos - 1
            end = pos
            varIndex = 1
        # Extract sequence and variant
        seq = extract_fasta(fasta, chrom, start, end)
        seqRef = seq[varIndex]
        # Check sequence equals reference
        if seqRef.upper() != ref:
            raise ValueError('reference {} variant {}'.format(name, seqRef))
        # Extract complexity and add to dataframe
        if any([x.islower() for x in seq]):
            if seqRef.islower():
                outputDF.loc[name] = 'within'
            else:
                outputDF.loc[name] = 'edge'
        else:
            outputDF.loc[name] = 'none'
    # Return data
    return(outputDF)

def homo_annotate(fasta, varList, flank = 100):
    # Create output dataframe
    varNames = [':'.join(map(str,x)) for x in varList]
    outputDF = pd.DataFrame(index = varNames, columns = ['monomer',
        'mono_ref', 'mono_var'])
    # Open fasta file and extract chromosome sizes
    fastaFile = pysam.FastaFile(fasta)
    chromSize = {}
    for chrom in fastaFile.references:
        chromSize[chrom] = fastaFile.get_reference_length(chrom)
    # Loop through indel list and extract annotation
    for name, (chrom, pos, ref, var) in zip(varNames, varList):
        # Skip non-indels
        if len(var) == 1 and '-' not in var:
            continue
        # Check variant is on chromosome
        if not 1 <= pos <= chromSize[chrom]:
            raise ValueError('position %s chromosome %s' %(pos, chrom))
        # Set start and end of interval
        start = max(1, pos - flank)
        end = min(chromSize[chrom], pos + flank)
        # Extract sequence
        sequence = extract_fasta(fasta, chrom, start, end)
        indelLocation = (pos - start)
        # Process indel
        if '-' in var:
            insertion = False
        else:
            insertion = True
        # Extract indel annotation and add to 
        indelData = indel_homopolymer(sequence, indelLocation, ref, var,
            insertion)
        outputDF.loc[name] = indelData
    # Return data
    return(outputDF)

def pairDist(read1, read2):
    ''' Function calculates inner and outer distance for paired reads if
    they are determined to be concordant. Function takes two arguemnts:
    
    1) read1 - Pysam AlignedSegment for read1.
    2) read2 - Pysam AlignedSegment for read2.
    
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

def chr_dict(bamFile):
    ''' Function returns a dictionary of chromosome data where the
    key is the target ID (tid) of the chromosome and the values are
    a tuple of the chromosome name and the chromosome length. Function
    takes the following 1 argument:
    
    1)  bamFile - Open Pysam AlignmentFile
    
    '''
    # Create tuple of chrosomes
    chrList = zip(bamFile.references, bamFile.lengths)
    # Add to dictionary
    chrDict = {}
    for c in chrList:
        tid = bamFile.gettid(c[0])
        chrDict[tid] = c
    # Return dictionary
    return(chrDict)
    
def startInterval(read, intervalSize, chrDict):
    ''' Function returns interval around the start of a read. If the
    read is aligned to forward strand then the start is defined as the
    end of the read aligned to the most 5' portion of the reference.
    If the read is aligned to the reverse strand then the start is defined
    as end of the read aligned to the most 3' portion of the reference.
    The function takes the following 3 arguments:
    
    1)  read - Pysam AlignedSegment
    2)  intervalSize - Size of the returned interval
    3)  chrDict - A dictionary where the keys are the target IDs (tid)
        of the chromosomes and the values is a tuple of the chromsome
        name and length.
    
    Function returns a tuple of the region in BED format; using zero
    based indexing and open ended interval and where name is read name
    and score is read mapping quality.
    
    '''
    # Extract flag
    flag = read.flag
    # Find strand of read and extract read start position
    if flag & 16:
        strand = '-'
        site = read.reference_end
    else:
        strand = '+'
        site = read.reference_start
    # Calculate start and end of interval
    halfSize = int(intervalSize / 2)
    start = site - halfSize
    end = site + halfSize
    # Extract chromosome data and adjust start and end
    chrom, length = chrDict[read.reference_id]
    start = max(0, start)
    end = min(length, end)
    # Extract mapping quality
    mapq = read.mapping_quality
    # Extract read name
    if flag & 64:
        name = read.query_name + '/1'
    elif flag & 128:
        name = read.query_name + '/2'
    else:
        name = read.query_name
    # Return data
    return((chrom, start, end, name, mapq, strand))

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
        read, groupdel = False
    ):
    # Extract cigar tuple and check
    cigartuple = read.cigartuples
    if cigartuple is None:
        return({})
    # Remove clipping from cigar
    cigartuple = [x for x in cigartuple if x[0] < 4]
    # Extract operations from tuple and check
    tupleSet = set([x[0] for x in cigartuple])
    if tupleSet.intersection([3, 6, 7, 8]):
        raise IOError('Irregular cigar found: %s' %(cigartuple))
    # Extract additional data
    sequence = list(read.query_alignment_sequence)
    quality = list(read.query_alignment_qualities)
    positions = read.get_reference_positions()
    # Process indels
    if tupleSet.intersection([1, 2]):
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
            # Group deletions
            if groupdel:
                # Alter sequence to include "-" to signify deletion
                sequence = sequence[:d[0]] + ['-' * d[1]] + sequence[d[0]:]
                # Alter quality to mean of flank
                dQual = (quality[d[0] - 1] + quality[d[0]]) / 2
                quality = quality[:d[0]] + [dQual] + quality[d[0]:]
                # Alter positions
                dPos = positions[d[0] - 1] + 1
                positions = positions[:d[0]] + [dPos] + positions[d[0]:]
            # Associate deletions with individual bases
            else:
                # Alter sequence to include "-" to signify deletion
                sequence = sequence[:d[0]] + ['-'] * d[1] + sequence[d[0]:]
                # Alter quality to mean of flank
                dQual = (quality[d[0] - 1] + quality[d[0]]) / 2
                quality = quality[:d[0]] + [dQual] * d[1] + quality[d[0]:]
                # Alter positions
                dPos = range(positions[d[0] - 1] + 1,
                    positions[d[0] - 1] + d[1] + 1)
                positions = positions[:d[0]] + dPos + positions[d[0]:]
    # Create and return output
    if len(sequence) == len(quality) and len(quality) == len(positions):
        positionDict = dict(zip(positions, zip(sequence, quality)))
        return(positionDict)
    else:
        raise IOError('Could not ascribe sequence to positions')

def extractPosition(
       openBam, chrom, position, minMapQ = 20, minBaseQ = 20, groupdel = False
    ):
    # check arguments
    toolbox.check_var(chrom, 'str')
    toolbox.check_var(position, 'int', mn = 1)
    toolbox.check_var(minMapQ, 'int', mn = 0)
    toolbox.check_var(minBaseQ, 'int', mn = 0)
    toolbox.check_var(groupdel, 'bool')
    # Set variables for mapping
    mapQuality = []
    baseCounts = {}
    # Loop through reads covering BAM
    for read in openBam.fetch(chrom, position-1, position):
        # Store mapping quality and skip reads with low values
        mapQuality.append(read.mapping_quality)
        if mapQuality[-1] < minMapQ:
            continue
        # Extract base calls for read
        baseDict = baseCalls(read, groupdel)
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

def extractPositionComplete(
       openBam, chrom, position, groupdel = False
    ):
    '''
    Function to extract all mapped base information from
    chromosomal position
    '''
    # check arguments
    toolbox.check_var(chrom, 'str')
    toolbox.check_var(position, 'int', mn = 1)
    toolbox.check_var(groupdel, 'bool')
    # Set variables for mapping
    baseCounts = {}
    # Loop through reads covering BAM
    for read in openBam.fetch(chrom, position-1, position):
        # Skip unmapped reads
        if read.is_unmapped:
            continue
        # Extract base calls for read
        baseDict = baseCalls(read, groupdel)
        # Extract reads for base of interest or skip
        try:
            base, baseQ = baseDict[position - 1]
        except KeyError:
            continue
        # Extract strand and mapping quality
        mapQ = read.mapping_quality
        strand = '-' if read.is_reverse else '+'
        # Add base to base dictionary
        if base in baseCounts:
            baseCounts[base].append((baseQ, mapQ, strand))
        else:
            baseCounts[base] = [(baseQ, mapQ, strand)]
    # Return data
    return(baseCounts)

def extractVariantFrequency(
        ref, var, baseData, minMapQ, minBaseQ
    ):
    # Create counter for data
    countDict = {ref : 0.0, var : 0.0}
    # Loop through base data
    for base, mapQ, baseQ in baseData:
        # Skip non reference or variant bases
        if base != ref or base != var:
            continue
        # Skip low mapping quality
        if mapQ < minMapQ:
            continue
        # Skip low base quality
        if baseQ < minBaseQ:
            continue
        # Add base to count
        countDict[base] += 1
    # Calculate and return frequency
    coverage = countDict[ref] + countDict[var]
    if coverage == 0:
        frequency = 0
    else:
        frequency = countDict[var] / coverage
    return(frequency)

def extractVariantCountsProcess(
        variantList, bamFile, pipe, minMapQ = 20, minBaseQ = 20,
        groupdel = False
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
    # check arguments
    toolbox.check_var(minMapQ, 'int', mn = 0)
    toolbox.check_var(minBaseQ, 'int', mn = 0)
    toolbox.check_var(groupdel, 'bool')
    # Create output dataframe
    variantNames = [':'.join(map(str, x)) for x in variantList]
    outData = pd.DataFrame(
        columns = ['refcount', 'reffor', 'varcount',
            'varfor', 'mapqual'],
        index = variantNames
    )
    # Open bamFile
    bam = pysam.AlignmentFile(bamFile)
    # Loop through variants and extract counts
    for name, (chrom, position, reference, variant) in zip(variantNames,
        variantList):
        # Extract base data
        baseCounts, mapqual = extractPosition(
            openBam = bam, chrom = chrom, position = position,
            minMapQ = minMapQ, minBaseQ = minBaseQ, groupdel = groupdel
        )
        # Extract reference counts
        if reference in baseCounts:
            refcount, reffor = baseCounts[reference]
        else:
            refcount = 0
            reffor = 0
        # Extract variant frequency
        if variant in baseCounts:
            varcount, varfor = baseCounts[variant]
        else:
            varcount = 0
            varfor = 0
        # Add data to output list
        outData.loc[name] = [refcount, reffor, varcount, varfor, mapqual]
    # Close bam
    bam.close()
    # Send data down pipe and close
    pipe.send(outData)
    pipe.close()

def calculateVariantMetrics(
        variantList, bamList, sampleNames, annovarPath, buildver, database,
        tempprefix, minMapQ = 20, minBaseQ = 20, groupdel = False,
        altQualNormal = None, homo = True, complexity = True, fasta = None
    ):
    ''' Function calculates metrics for variants across multiple samples

    Args:
        variantList (list): A list of four element tuples that list the chrom
            position, reference and variant.
        bamList (list): A list of BAM files from which to extract variant
            annotation.
        sampleNames (list): A list of sample names for each of the BAM files.
        annovarPath (str): Path to annovar executable
        buildver (str): Genome build to use for annotation.
        database (str): Database to use for annotation.
        homo (bool): Whether to annotate indels for overlapping homopolymers.
        complexity (bool): Whether to annotate variants for complexity using
            soft-masking in FASTA file.
        fasta (str): Full path to FASTA file. Required for  hompolymer
            annotation.
    
    '''
    # Check arguments
    toolbox.check_var(minMapQ, 'int', mn = 0)
    toolbox.check_var(minBaseQ, 'int', mn = 0)
    toolbox.check_var(groupdel, 'bool')
    toolbox.check_var(altQualNormal, 'int', mn = 2)
    # Check supplied Names
    if not isinstance(sampleNames, (list, tuple)):
        raise IOError('sampleNasmes must be a list or a tuple')
    if len(sampleNames) != len(bamList):
        raise IOError('Must be a sample name for each BAM')
    # Create output dataframe
    varnames = [':'.join(map(str,x)) for x in variantList]
    outputData = pd.DataFrame(index = varnames, columns = ['chr','pos','ref',
        'var','minp'])
    for x in varnames:
        chrom, position, ref, var = x.split(':')
        outputData.loc[x] = [chrom, int(position), ref, var, 1]
    # Add homopolymer annotation
    if homo:
        homoData = homo_annotate(fasta, variantList, flank = 100)
        outputData = pd.concat([outputData, homoData], axis = 1)
    if complexity:
        compData = comp_annotate(fasta, variantList)
        outputData = pd.concat([outputData, compData], axis = 1)
    # Create variables to store pipe and process data
    processDict = {}
    # Create process for each BAM
    for number, (name, bamFile) in enumerate(zip(sampleNames, bamList)):
        # Process reference sample
        if number == 0:
            # Create pipes and process
            pipeRecv, pipeSend = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = extractVariantCountsProcess,
                args = (variantList, bamFile, pipeSend, minMapQ, minBaseQ,
                    groupdel)
            )
            process.start()
            pipeSend.close()
            # Store data
            processDict['reference'] = (name, pipeRecv, process)
            # Add extra process for reduced base quality in normal
            if altQualNormal:
                # Create pipes and process
                pipeRecv, pipeSend = multiprocessing.Pipe(False)
                process = multiprocessing.Process(
                    target = extractVariantCountsProcess,
                    args = (variantList, bamFile, pipeSend, minMapQ,
                        altQualNormal, groupdel)
                )
                process.start()
                pipeSend.close()
                # Store data
                processDict['altref'] = (name, pipeRecv, process)
        # Process non-reference samples
        else:
            # Create pipes and process
            pipeRecv, pipeSend = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = extractVariantCountsProcess,
                args = (variantList, bamFile, pipeSend, minMapQ, minBaseQ,
                    groupdel)
            )
            process.start()
            pipeSend.close()
            # Store data
            processDict[number] = (name, pipeRecv, process)
    # Extract data for reference
    name, pipe, process = processDict.pop('reference')
    sampleData = pipe.recv()
    pipe.close()
    process.join()
    # Add data for reference to output
    outputData[name + '_ref'] = sampleData['refcount']
    outputData[name + '_var'] = sampleData['varcount']
    outputData[name + '_freq'] = sampleData['varcount'] / (
        sampleData['refcount'] + sampleData['varcount'])
    outputData[name + '_mapq'] = sampleData['mapqual']
    # Store data for reference
    normRef = outputData[name + '_ref']
    normVar = outputData[name + '_var']
    normFreq = outputData[name + '_freq']
    # Extract data for alternative frequence
    if 'altref' in processDict:
        # Extract data for alternative frequency
        name, pipe, process = processDict.pop('altref')
        sampleData = pipe.recv()
        pipe.close()
        process.join()
        # Add data for alternative frequency to output
        outputData[name + '_altfreq'] = sampleData['varcount'] / (
            sampleData['refcount'] + sampleData['varcount'])
    # Extract variants for each BAM
    for key in processDict:
        # Extract data for reference
        name, pipe, process = processDict[key]
        sampleData = pipe.recv()
        pipe.close()
        process.join()
        # Add data to output
        outputData[name + '_ref'] = sampleData['refcount']
        outputData[name + '_var'] = sampleData['varcount']
        outputData[name + '_freq'] = sampleData['varcount'] / (
            sampleData['refcount'] + sampleData['varcount'])
        outputData[name + '_mapq'] = sampleData['mapqual']
        # Calculate pvalue
        pvalue = []
        for freq, normal, sample in zip(
            zip(normFreq, outputData[name + '_freq']),
            zip(normRef, normVar), 
            zip(outputData[name + '_ref'], outputData[name + '_var'])
        ):
            if freq[0] < freq[1]:
                pvalue.append(fisher_exact(
                    [normal, sample],
                    alternative = 'greater'
                )[1])
            elif freq[0] > freq[1]:
                pvalue.append(fisher_exact(
                    [normal, sample],
                    alternative = 'less'
                )[1])
            else:
                pvalue.append(fisher_exact(
                    [normal, sample],
                    alternative = 'two-sided'
                )[1])
        # Rename tables columns and append to output
        outputData[name + '_pvalue'] = pvalue
    # Calculate minium pvalue and sort
    pvalueIndex = [x.endswith('pvalue') for x in outputData.columns]
    minp = outputData.loc[:,pvalueIndex].min(1)
    outputData['minp'] = minp
    # Add annovar annotation and concat to dataframe
    geneAnno = annovar.geneAnno2DF(variantList = variantList,
        path = annovarPath, buildver = buildver, database = database,
        tempprefix = tempprefix)
    outputData = pd.concat([outputData, geneAnno], axis = 1)
    outputData.sort_values('minp', inplace = True)
    # Return data
    return(outputData)

def filterVariantMetrics(
        inFile, outFile, minNormCov = 8, minTumrCov = 10, maxNormFreq = 0.01,
        minTumrFreq = 0.1, maxPvalue = 0.01, minTumrVar = 2, minMapQ = 40,
        rmIndels = False
    ):
    # Generate counter for filtering
    logData = collections.OrderedDict([
        ('Total', 0),
        ('Indels', 0),
        ('Normal coverage', 0),
        ('Normal frequency', 0),
        ('Normal mapping', 0),
        ('Tumour coverage', 0),
        ('Tumour frequency', 0),
        ('Tumour mapping', 0),
        ('Tumour pvalue', 0),
        ('Accepted', 0)
    ])
    # Open input file and extract header
    inputf = open(inFile)
    header = inputf.next().strip().split('\t')
    # Check expected header start and end
    if tuple(header[0:5]) != ('chr', 'pos', 'ref', 'var', 'minp'):
        raise ValueError('Unexpected header')
    if tuple(header[-4:]) != ('class', 'genes', 'affect', 'change'):
        raise ValueError('Unexpected header')
    # Extract indices for normal
    if tuple(header[5:9]) == ('monomer', 'mono_ref', 'mono_var', 'low_comp'):
        normInd = slice(9,13)
    elif tuple(header[5:8]) == ('monomer', 'mono_ref', 'mono_var'):
        normInd = slice(8,12)
    elif header[4] == 'low_comp':
        normInd = slice(6,10)
    else:
        normInd = slice(5,9)
    # Extract indices for tumour
    tumrStart = normInd.stop
    if header[tumrStart].endswith('altfreq'):
        tumrStart += 1
    tumrIndList = []
    for x in range(tumrStart, len(header) - 4, 5):
        tumrIndList.append(slice(x, x + 5))
    # Check indices for normal and tumours
    normRef, normVar, normFreq, normMapq = header[normInd]
    if (not normRef.endswith('_ref') or
        not normVar.endswith('_var') or
        not normFreq.endswith('_freq') or
        not normMapq.endswith('_mapq')):
        raise ValueError('Unexpected header')
    for tumrInd in tumrIndList:
        tumrRef, tumrVar, tumrFreq, tumrMapq, tumrPv = header[tumrInd]
        if (not tumrRef.endswith('_ref') or
            not tumrVar.endswith('_var') or
            not tumrFreq.endswith('_freq') or
            not tumrMapq.endswith('_mapq') or
            not tumrPv.endswith('_pvalue')):
            raise ValueError('Unexpected header')
    # Calculate tumour number
    tumourNo = float(len(tumrIndList))
    # Open output file and add header
    outputf = open(outFile, 'w')
    outputf.write('\t'.join(header) + '\n')
    # Loop through input file
    for line in inputf:
        # Count and split line
        logData['Total'] += 1
        lineData = line.strip().split('\t')
        # Perform indel filtering
        if rmIndels:
            varSeq = lineData[3]
            if '-' in varSeq or len(varSeq) > 1:
                logData['Indels'] += 1
                continue
        # Extract and check normal
        normRef, normVar, normFreq, normMapq = lineData[normInd]
        normCov = int(normRef) + int(normVar)
        if normCov < minNormCov:
            logData['Normal coverage'] += 1
            continue
        if float(normFreq) > maxNormFreq:
            logData['Normal frequency'] += 1
            continue
        if int(normMapq) < minMapQ:
            logData['Normal mapping'] += 1
            continue
        # Identify acceptable tumours
        tumrPass = False
        tumrLog = collections.defaultdict(int)
        tumrPvalue = []
        for tumrInd in tumrIndList:
            # Extract tumour data
            tumrRef, tumrVar, tumrFreq, tumrMapq, tumrPv = lineData[tumrInd]
            tumrCov = int(tumrRef) + int(tumrVar)
            # Count and skip tumour if tumour data unacceptable
            if tumrCov < minTumrCov:
                tumrLog['Tumour coverage'] += 1
                lineData[tumrInd.stop - 1] = ''
                continue
            if float(tumrFreq) < minTumrFreq:
                tumrLog['Tumour frequency'] += 1
                lineData[tumrInd.stop - 1] = ''
                continue
            if int(tumrMapq) < minMapQ:
                tumrLog['Tumour mapping'] += 1
                lineData[tumrInd.stop - 1] = ''
                continue
            if float(tumrPv) > maxPvalue:
                tumrLog['Tumour pvalue'] += 1
                lineData[tumrInd.stop - 1] = ''
                continue
            # Process data if tumour acceptable
            tumrPvalue.append(float(tumrPv))
            tumrPass = True
        # Count and extract accepted variants
        if tumrPass:
            logData['Accepted'] += 1
            lineData[4] = str(min(tumrPvalue))
            outputf.write('\t'.join(lineData) + '\n')
        # Else log reason for tumour dismissal
        else:
            for key in tumrLog:
                logData[key] += (tumrLog[key] / tumourNo)
        # Count and output accepted variants
    # Close files and return log
    inputf.close()
    outputf.close()
    return(logData)

def strandTranAlign(alignment):
    ''' Function calculates the proportion of specificity of reads aligned to a
    sense reference transcriptome. Function takes 1 arguments:
    
    1)  alignment - A sam/bam file. Sorting unrequired.
    
    '''
    # Open alignment file and create variable to store count
    alignFile = pysam.AlignmentFile(alignment)
    count = {'total': 0, 'sense': 0, 'antisense': 0, 'discord': 0}
    # Loop through reads in file
    for read in alignFile:
        # Extract flag
        flag = read.flag
        # Skip unmapped reads
        if flag & 4:
            continue
        # Count mapped reads
        else:
            count['total'] += 1
        # Skip improper pairs
        if not flag & 2:
            count['discord'] += 1
            continue
        # Process read1
        if flag & 64:
            # Process R1 reverse
            if flag & 16:
                # Process R1 + R2 reverse
                if flag & 32:
                    count['discord'] += 1
                # Process R1 reverse R2 forward
                else:
                    count['antisense'] += 1
            # Process R1 forward
            else:
                # Process R1 forward + R2 reverse
                if flag & 32:
                    count['sense'] += 1
                # Process R1 + R2 forward
                else:
                    count['antisense'] += 1
        # Process R2
        elif flag & 128:
            # Process R2 reverse
            if flag & 16:
                # Process R2 + R1 reverse
                if flag & 32:
                    count['discord'] += 1
                # Process R2 reverse R1 forward
                else:
                    count['sense'] += 1
            # Process R2 forward
            else:
                # Process R2 forward + R1 reverse
                if flag & 32:
                    count['antisense'] += 1
                # Process R2 + R1 forward
                else:
                    count['discord'] += 1
    # Retrun counts
    return(counts)

