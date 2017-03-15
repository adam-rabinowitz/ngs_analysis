from general_python import toolbox
import collections
import bisect

def calcRatio(
        mpileup1, mpileup2
    ):
    # Create an return awk command
    awkCount = 'awk \'BEGIN{FS = "\\t"};{SUM += $4};END{print SUM}\''
    finalCommand = 'awk -v N=$(%s %s) -v T=$(%s %s) \'BEGIN{print N / T}\'' %(
        awkCount, mpileup1, awkCount, mpileup2
    )
    return(finalCommand)

def germlineSNP(
        mpileup, outfile, minCov = 8, minVarReads = 2, minAvgQual = 15,
        minVarFreq = 0.01, minHomFreq = 0.75, pValue = 0.99, javaPath = 'java',
        varscanPath = 'varscan.jar', memory = None
    ):
    ''' A function to generate a varscan command to identify germline snps.
    
    Args:
        mpileup (str)- Full path to mpileup file.
        outfile (str)- Full path to output file.
        minCov (int)- Min coverage to call variant.
        minVarReads (int)- Min number of variant reads to call variant.
        minAvgQual (int)- Min average base qaulity to call variant.
        minVarFreq (float/int)- Min variant frequency.
        pValue (float)- Minimum pvalue to call vatiant
        javaPath (str)- Java executable.
        varscanPath (str)- Path to java jar file.
        memory (int)- Memory in GB for java stack.
    
    '''
    # Check arguments
    if not isinstance(minCov, int):
        raise TypeError('minCov must be integer')
    if 1 > minCov:
        raise ValueError('minCov must be >0')
    if not isinstance(minVarReads, int):
        raise TypeError('minVarReads must be integer')
    if 1 > minVarReads:
        raise ValueError('minVarReads must be >0')
    if not isinstance(minAvgQual, int):
        raise TypeError('minAvgQual must be integer')
    if 1 > minAvgQual:
        raise ValueError('minVarReads must be >0')
    if not isinstance(minVarFreq, (float, int)):
        raise TypeError('minVarFreq must be float or integer')
    if not 0 < minVarFreq <= 1:
        raise ValueError('minVarFreq must be >0 and <=1')
    if not isinstance(minHomFreq, (float, int)):
        raise TypeError('minHomFreq must be float or integer')
    if not 0 < minHomFreq <=1:
        raise ValueError('minHomFreq must be >0 and <=1')
    if not isinstance(pValue, float):
        raise TypeError('pValue must be float or integer')
    if not 0 < pValue <1:
        raise ValueError('pValue must be >0 and <1')
    if not isinstance(javaPath, str):
        raise TypeError('javaPath must be string')
    if not isinstance(varscanPath, str):
        raise TypeError('pValue must be string')
    if memory and not isinstance(memory, int):
        raise TypeError('memory must be a string')
    if memory and not 0 < memory <=256:
        raise ValueError('memory must be >0 and <=256')
    # Create command
    command = [
        javaPath, '-jar', varscanPath, 'pileup2snp', mpileup, '--min-coverage',
        str(minCov), '--min-reads2', str(minVarReads), '--min-avg-qual',
        str(minAvgQual), '--min-var-frequency', str(minVarFreq),
        '--min-freq-for-hom', str(minHomFreq), '--p-values', str(pValue),
        '>', outfile
    ]
    # Add memory
    if memory:
        command.insert(1, '-Xmx{}g'.format(memory))
    # Join and return command
    command = ' '.join(command)
    return(command)

def filterGermlineSNP(
        infile, outfile, minCov = 10, minVarReads = 2, minVarFreq = 0.1,
        maxVarFreq = 1, maxStrandBias=0.9
    ):
    # Check arguments
    if not isinstance(minCov, int):
        raise TypeError('minCov must be integer')
    if 1 > minCov:
        raise ValueError('minCov must be >0')
    if not isinstance(minVarReads, int):
        raise TypeError('minVarReads must be integer')
    if 1 > minVarReads:
        raise ValueError('minVarReads must be >0')
    if not isinstance(minVarFreq, (float, int)):
        raise TypeError('minVarFreq must be float or integer')
    if not 0 < minVarFreq <= 1:
        raise ValueError('minVarFreq must be >0 and <=1')
    if not isinstance(maxVarFreq, (float, int)):
        raise TypeError('maxVarFreq must be float or integer')
    if not minVarFreq <= maxVarFreq <= 1:
        raise ValueError('maxVarFreq must be >minVarFreq and <=1')
    # Open input and output files
    ifile = open(infile, 'r')
    ofile = open(outfile, 'w')
    # Extract and process header
    header = ifile.next()
    ofile.write(header)
    header = header.strip().split('\t')
    refReadsIndex = header.index('Reads1')
    varReadsIndex = header.index('Reads2')
    varForStrandIndex = header.index('Reads2Plus')
    # Loop though input file
    for line in ifile:
        # Extract data
        linedata = line.strip().split('\t')
        refReads = float(linedata[refReadsIndex])
        varReads = float(linedata[varReadsIndex])
        forStrand = float(linedata[varForStrandIndex])
        # Perform calculations
        cov = refReads + varReads
        varFreq = varReads / (refReads + varReads)
        strandBias = forStrand / varReads
        strandBias = max(1 - strandBias, strandBias)
        # Filter
        if minCov > cov:
            continue
        if minVarReads > varReads:
            continue
        if minVarFreq > varFreq:
            continue
        if maxVarFreq < varFreq:
            continue
        if strandBias > maxStrandBias:
            continue
        # Write output
        ofile.write(line)
    # Close files
    ifile.close()
    ofile.close()

def somatic(
        mpileup1, mpileup2, outPrefix, purity = 0.5, minCovNormal = 8,
        minCovTumour = 6, minHetFreq = 0.1, minHomFreq = 0.75,
        normalPurity = 1.0, tumourPurity = 0.5, pValueHet = 0.99,
        pValueSomatic = 0.05, strandFilter = False, javaPath = 'java',
        varscanPath = 'varscan.jar', memory = None
    ):
    # Check commands
    if not isinstance(purity, (float, int)):
        raise TypeError('purity must be float or integer')
    if not 0 < purity <= 1:
        raise ValueError('purity must be >0 and <=1')
    if not isinstance(minCovNormal, int):
        raise TypeError('minCovNormal must be an integer')
    if minCovNormal <= 0:
        raise ValueError('minCovNormal must be >0')
    if not isinstance(minCovTumour, int):
        raise TypeError('minCovTumour must be an integer')
    if minCovTumour <= 0:
        raise ValueError('minCovTumour must be >0')
    if not isinstance(minHetFreq, float):
        raise TypeError('minHetFreq must be float')
    if not 0 < minHetFreq < 1:
        raise ValueError('minHetFreq must be >0 and <1')
    if not isinstance(minHomFreq, (float,int)):
        raise TypeError('minHomFreq must be float')
    if not 0 < minHomFreq <= 1:
        raise ValueError('minHomFreq must be >0 and <=1')
    if not isinstance(normalPurity, (float, int)):
        raise TypeError('normalPurity must be float or integer')
    if not 0 < normalPurity <= 1:
        raise ValueError('purity must be >0 and <=1')
    if not isinstance(pValueHet, (float, int)):
        raise TypeError('pValueHet must be float or integer')
    if not 0 <= pValueHet <= 1:
        raise ValueError('pValueHet must be >=0 and <=1')
    if not isinstance(pValueSomatic, (float, int)):
        raise TypeError('pValueSomatic must be float or integer')
    if not 0 <= pValueSomatic <= 1:
        raise ValueError('pValueSomatic must be >=0 and <=1')
    if not isinstance(strandFilter, bool):
        raise TypeError('strandFilter must be bool')
    if memory and not isinstance(memory, int):
        raise TypeError('memory must be an integer')
    if memory and not 1 <= memory <= 256:
        raise ValueError('memory must be >=1 and <=256')
    # Create command
    command = [
        javaPath, '-jar', varscanPath, 'somatic', mpileup1, mpileup2,
        outPrefix, '--min-coverage-normal', str(minCovNormal),
        '--min-coverage-tumor', str(minCovTumour), '--min-var-freq',
        str(minHetFreq), '--min-freq-for-hom', str(minHomFreq),
        '--normal-purity', str(normalPurity), '--tumor-purity',
        str(tumourPurity), '--p-value', str(pValueHet), '--somatic-p-value',
        str(pValueSomatic)
    ]
    if strandFilter:
        command.extend(['--strand-filter', '1'])
    # Add memory
    if memory:
        command.insert(1, '-Xmx{}g'.format(memory))
    # Join and return command
    command = ' '.join(command)
    return(command)

def filterVarscan(
        inFile, outFile, filterFile = None,  minCovNormal = 10,
        minCovTumour = 10, minFreqTumour = 0.05, maxFreqNormal = 1,
        minVarTumour = 2, maxPvalue = 0.05, somatic = True, flank = 25,
        maxNeighbour = 0
    ):
    # Create counter
    logData =collections.OrderedDict([
        ('Total', 0),
        ('Somatic status', 0),
        ('P-value', 0),
        ('Tumour coverage', 0),
        ('Tumour frequency', 0),
        ('Tumour count', 0),
        ('Normal coverage', 0),
        ('Normal frequency', 0),
        ('Neighbours', 0),
        ('Passed filters', 0)
    ])
    # Check variables
    toolbox.check_var(inFile, 'file')
    toolbox.check_var(filterFile, 'file')
    toolbox.check_var(minCovNormal, 'int', mn = 1)
    toolbox.check_var(minCovTumour, 'int', mn = 1)
    toolbox.check_var(minFreqTumour, 'num', gt = 0, mx = 1)
    toolbox.check_var(maxFreqNormal, 'num', mn = 0, mx = 1)
    toolbox.check_var(minVarTumour, 'int', mn = 1)
    toolbox.check_var(maxPvalue, 'num', gt = 0)
    toolbox.check_var(somatic, 'bool')
    toolbox.check_var(flank, 'int', mn = 0)
    toolbox.check_var(maxNeighbour, 'int', mn = 0)
    # Create dictionary to store variant positions
    varPos = {}
    # Extract coordinates for neighbour filtering
    for varFile in [inFile, filterFile]:
        if varFile is None:
            continue
        with open(varFile) as varIn:
            header = varIn.next()
            for line in varIn:
                chrom, pos = line.split('\t')[:2]
                if chrom in varPos:
                    varPos[chrom].append(int(pos))
                else:
                    varPos[chrom] = [int(pos)]
    # Sort data
    for key in varPos:
        varPos[key].sort()
    # Open input and output files
    with open(inFile) as varin:
        with open(outFile, 'w') as varout:
            # Write header
            varout.write(varin.next())
            # Loop through input
            for line in varin:
                # Count and extract data
                logData['Total'] += 1
                varData = line.split('\t')
                # Check somatic status and p-value
                status = str(varData[12])
                pValue = float(varData[14])
                if somatic and status != 'Somatic':
                    logData['Somatic status'] += 1
                    continue
                if pValue > maxPvalue:
                    logData['P-value'] += 1
                    continue
                # Check coverage and frequency
                covNormal = int(varData[4]) + int(varData[5])
                freqNormal = int(varData[5]) / float(covNormal)
                covTumour = int(varData[8]) + int(varData[9])
                freqTumour = int(varData[9]) / float(covTumour)
                varTumour = int(varData[9])
                if covTumour < minCovTumour:
                    logData['Tumour coverage'] += 1
                    continue
                if freqTumour < minFreqTumour:
                    logData['Tumour frequency'] += 1
                    continue
                if varTumour < minVarTumour:
                    logData['Tumour count'] += 1
                    continue
                if covNormal < minCovNormal:
                    logData['Normal coverage'] += 1
                    continue
                if freqNormal > maxFreqNormal:
                    logData['Normal frequency'] += 1
                    continue
                # Check flanking mutations
                chrom = varData[0]
                start = int(varData[1]) - flank
                end = int(varData[1]) + flank
                startIndex = bisect.bisect_left(varPos[chrom], start)
                endIndex = bisect.bisect_right(varPos[chrom], end, lo = startIndex)
                neighbourCount = (endIndex - startIndex) - 1
                if neighbourCount > maxNeighbour:
                    logData['Neighbours'] += 1
                    continue
                # Write output line
                logData['Passed filters'] += 1
                varout.write(line)
    # Return log
    return(logData)

def copynumber(
        mpileup1, mpileup2, outPrefix, minBaseQ = 20, minMapQ = 20,
        minCov = 20, minSegSize = 10, maxSegSize = 100, pValue = 0.01,
        dataRatio = None, javaPath = 'java', varscanPath = 'varscan.jar',
        memory = None
    ):
    ''' Function to generate command to perform copynumber calling using
    the varscan program. Function takes the following X arguments:
    
    1)  pileup - The normal-tumour pileup.
    2)  outPrefix - The prefix of the output files.
    3)  minBaseQ - Minimum base quality for coverage.
    4)  minMapQ - Minimum read mapping quality for coverage.
    5)  minCov - Minimum coverage for copynumber segments.
    6)  minSegSize - Minimum segment size.
    7)  maxSegSize - Maximum segment size.
    8)  pValue - P-value for significant copynumber change-point.
    9)  dataRatio - The normal/tumor input data ratio.
    
    '''
    # Check commands
    if not isinstance(minBaseQ, int):
        raise TypeError('minBaseQ must be integer')
    if minBaseQ < 1:
        raise ValueError('minBaseQ must be >0')
    if not isinstance(minMapQ, int):
        raise TypeError('minMapQ must be integer')
    if minMapQ < 1:
        raise ValueError('minMapQ must be >0')
    if not isinstance(minCov, int):
        raise TypeError('minCov must be integer')
    if minCov < 1:
        raise ValueError('minCov must be >0')
    if not isinstance(minSegSize, int):
        raise TypeError('minSegSize must be integer')
    if minSegSize < 1:
        raise ValueError('minSegSize must be >0')
    if not isinstance(maxSegSize, int):
        raise TypeError('maxSegSize must be integer')
    if maxSegSize < minSegSize:
        raise ValueError('maxSegSize must be >=minSegSize')
    if not isinstance(pValue, (int, float)):
        raise TypeError('pValue must be integer or float')
    if not 0 <= pValue <= 1:
        raise ValueError('pValue must be >=0 and <=1')
    if dataRatio and not isinstance(dataRation (int, float)):
        raise TypeError('dataRatio must be float or integer')
    if dataRatio and not 0.01 <= dataRatio <= 100:
        raise ValueError('dataRatio must be >=0.01 and <=100')
    if memory and not isinstance(memory, int):
        raise TypeError('memory must be an integer')
    if memory and not 1 <= memory <= 256:
        raise ValueError('memory must be >=1 and <=256')
    # Create command to calculate depth if required
    if dataRatio is None:
        ratioCommand = 'R=$(%s) && echo "Ratio: $R"' %(
            calcRatio(mpileup1, mpileup2))
        dataRatio = '$R'
    else:
        ratioCommand = ''
    # Create copy number command
    copyCommand = [
        javaPath, '-jar', varscanPath, 'copynumber', mpileup1, mpileup2,
        outPrefix, '--min-base-qual', str(minBaseQ), '--min-map-qual',
        str(minMapQ), '--min-coverage', str(minCov), '--min-segment-size',
        str(minSegSize), '--max-segment-size', str(maxSegSize), '--p-value',
        str(pValue), '--data-ratio', str(dataRatio)
    ]
    # Add memory
    if memory:
        copyCommand.insert(1, '-Xmx{}g'.format(memory))
    # Combine and return commands
    if ratioCommand:
        return('%s && %s' %(ratioCommand, ' '.join(copyCommand)))
    else:
        return(' '.join(copyCommand))

def filterSomatic(
        inFile, outFile, minCov = 10, minReads = 2, minStrands = 1,
        minAvgQ = 10, minVarFreq = 0.1, pValue = 0.05, indelFile = None,
        javaPath = 'java', varscanPath = 'varscan.jar'
    ):
    '''
    1)  minCov - Minimum read depth.
    2)  minReads - Minimum supporting reads for a variant.
    3)  minStrands - Minimum number of strands on which variant observed.
    4)  minAvgQ - Minimum average base quality for variant-supporting reads.
    5)  minVarFreq - Minimum variant allele frequency threshold.
    6)  pValue - Default p-value threshold for calling variants.
    7)  indelFile - File of indels for filtering nearby SNPs.
    8)  outFile - Output file for filtered variants.
    
    '''
    # Check numerical arguments
    toolbox.check_var(minCov, 'int', mn = 1)
    toolbox.check_var(minReads, 'int', mn = 1)
    toolbox.check_var(minStrands, 'int', mn = 1, mx = 2)
    toolbox.check_var(minAvgQ, 'int', mn = 2)
    toolbox.check_var(minVarFreq, 'num', gt = 0, mx = 1)
    toolbox.check_var(pValue, 'num', gt = 0, mx = 1)
    # Create command
    command = [javaPath, '-jar', varscanPath, 'somaticFilter', inFile,
        '--min-coverage', str(minCov), '--min-reads2', str(minReads),
        '--min-strands2', str(minStrands), '--min-avg-qual', str(minAvgQ),
        '--min-var-freq', str(minVarFreq), '--p-value', str(pValue),
        '--output-file', outFile]
    # Append indel file if supplied
    if indelFile:
        command.extend(['--indel-file', indelFile])
    # Return command
    command = ' '.join(command)
    return command

def variantSet(
        varscanFile, somatic = False
    ):
    '''
    Function to convert an varscan putput file to an annovar input file
    
    '''
    # Create output list
    varSet = set()
    # Open input file and skip header
    with open(varscanFile, 'r') as inputFile:
        # Process header
        inputHeader = inputFile.next().strip().split('\t')
        somaticIndex = inputHeader.index('somatic_status')
        # Loop through file and save data as av input
        for line in inputFile:
            # extract line data
            lineData = line.strip().split('\t')
            # Handle somatic status
            if somatic and lineData[somaticIndex] != 'Somatic':
                continue
            # Extract data
            chrom, position, reference, variants = lineData[:4]
            position = int(position)
            # Loop through potential variants
            for var in variants.split('/'):
                # Process SNV
                if len(var) == 1:
                    # Store variant data
                    varSet.add((chrom, position, reference, var))
                # Process deletions
                elif "-" in var:
                    # Modify positional, reference and variant values
                    adjpos = position + 1
                    adjref = var[1]
                    adjvar = "-" * (len(var) - 1)
                    # Store variant data
                    varSet.add((chrom, adjpos, adjref, adjvar))
                # Process insertions
                elif '+' in var:
                    # Modify variant values
                    adjvar = reference + var[1:]
                    # Store variant data
                    varSet.add((chrom, position, reference, adjvar))
                # Raise error for unrecognised variant
                else:
                    raise IOError('Variant %s could not be parsed' %(var))
    # Close input and return variant list
    return(varSet)

#varscan2annovar(
#    inFile = '/farm/scratch/rs-bio-lif/rabino01/Elza/varscan/310123-T1R1.filter.somatic.snp',
#    outFile = '/farm/scratch/rs-bio-lif/rabino01/Elza/310123-T1R1.filter.somatic.av'
#)
