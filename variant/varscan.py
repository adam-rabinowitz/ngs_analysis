from general_python import toolbox

def calcRatio(
        mpileup1, mpileup2
    ):
    # Create an return awk command
    awkCount = 'awk \'BEGIN{FS = "\\t"};{SUM += $4};END{print SUM}\''
    finalCommand = 'awk -v N=$(%s %s) -v T=$(%s %s) \'BEGIN{print N / T}\'' %(
        awkCount, mpileup1, awkCount, mpileup2
    )
    return(finalCommand)

def somatic(
        mpileup1, mpileup2, outPrefix, purity = 0.5, minCovNormal = 8,
        minCovTumour = 6, minHetFreq = 0.1, minHomFreq = 0.75,
        normalPurity = 1.0, tumourPurity = 0.5, pValueHet = 0.99,
        pValueSomatic = 0.05, strandFilter = False, javaPath = 'java',
        varscanPath = 'varscan.jar'
    ):
    # Check commands
    toolbox.checkArg(purity, 'num', gt = 0, mx = 1)
    toolbox.checkArg(minCovNormal, 'int', gt = 0)
    toolbox.checkArg(minCovTumour, 'int', gt = 0)
    toolbox.checkArg(minHetFreq, 'num', gt = 0, lt = 1)
    toolbox.checkArg(minHomFreq, 'num', gt = 0, mx = 1)
    toolbox.checkArg(normalPurity, 'num', gt = 0, mx = 1)
    toolbox.checkArg(normalPurity, 'num', gt = 0, mx = 1)
    toolbox.checkArg(pValueHet, 'num', mn = 0, mx = 1)
    toolbox.checkArg(pValueSomatic, 'num', mn = 0, mx = 1)
    toolbox.checkArg(strandFilter, 'bool')
    toolbox.checkArg(javaPath, 'file')
    toolbox.checkArg(varscanPath, 'file')
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
    # Join and return command
    command = ' '.join(command)
    return(command)

def copynumber(
        mpileup1, mpileup2, outPrefix, minBaseQ = 20, minMapQ = 20,
        minCov = 20, minSegSize = 10, maxSegSize = 100, pValue = 0.01,
        dataRatio = None, javaPath = 'java', varscanPath = 'varscan.jar'
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
    toolbox.checkArg(minBaseQ, 'int', gt = 0)
    toolbox.checkArg(minMapQ, 'int', gt = 0)
    toolbox.checkArg(minCov, 'int', gt = 0)
    toolbox.checkArg(minSegSize, 'int', gt = 0)
    toolbox.checkArg(maxSegSize, 'int', mn = minSegSize)
    toolbox.checkArg(pValue, 'num', mn = 0, mx = 1)
    toolbox.checkArg(dataRatio, 'num', mn = 0.01, mx = 100)
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
    # Combine and return commands
    if ratioCommand:
        return('%s && %s' %(ratioCommand, ' '.join(copyCommand)))
    else:
        return(' '.join(copyCommand))
