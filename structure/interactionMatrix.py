import collections
import multiprocessing
import os
import tempfile
import subprocess
import gzip
import cStringIO
import numpy as np
from ngs_python.structure import analyseInteraction

class genomeBin(object):

    def __init__(self, binData):
        self.binData = binData
        if isinstance(self.binData, str):
            self.binDict, self.binNames, self.binCount = self.readBed()
        elif isinstance(self.binData, tuple):
            self.binDict, self.binNames, self.binCount = self.createBins()
        else:
            raise IOError('tuple or string must be supplied')
    
    def createBins(self):
        # Unpack tuple
        chrFile, maxWidth, fixedWidth = self.binData
        # Create output data variables
        binDict = collections.OrderedDict()
        binNames = []
        count = 0
        # Extract chromosome name and length
        chrData = []
        with open(chrFile) as inFile:
            for line in inFile:
                chrName, chrLength = line.strip().split('\t')
                chrData.append((chrName, chrLength))
        # Loop through chromsome and generate bin data frame
        for chrName, chrLength in chrData:
            # Generate bin start and end for equal width bins
            if fixedWidth:
                # Calculate bin number
                remainder = float(int(chrLength) % maxWidth)
                # Calculate start and stop of first bin
                start = int(np.floor(remainder/2)) 
                end = start + maxWidth
                # Calculate start and stop of remianing bins
                binEnd = np.arange(end, int(chrLength) + 1, maxWidth,
                    dtype = np.uint32)
                binStart = (binEnd - maxWidth) + 1
            # Generate bin start and end sites for unequal width bins
            else:
                # Calculate number of bins
                binNo = np.ceil(float(chrLength) / maxWidth)
                # Calculate excess number of elements within bins
                excess = (binNo * maxWidth) - int(chrLength)
                # Calculate bin size 
                largeBin = int(maxWidth - np.floor(excess / binNo))
                smallBin = int(maxWidth - np.ceil(excess / binNo))
                # Calculate bin numbers
                smallBinNo = int(excess % binNo )
                largeBinNo = int(binNo - smallBinNo)
                # Generate bin width
                binWidth = np.array([largeBin] * largeBinNo + 
                    [smallBin] * smallBinNo)
                # Generate bins
                binEnd = np.cumsum(binWidth, dtype = np.uint32)
                binStart = (binEnd - binWidth) + 1
            # Generate arrays containing bin index and chromosome name
            binIndex = np.arange(count, count + len(binStart),
                dtype = np.uint32)
            count += len(binStart)
            binChr = np.array([chrName] *  len(binStart))
            # Store output data
            for bin in np.transpose([binChr, binStart, binEnd]):
                binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
            binDict[chrName] = {'start' : binStart, 'end' : binEnd,
                'index' : binIndex, 'count' : len(binStart)}
        # Return dataframe
        return(binDict, binNames, count)
    
    def readBed(self):
        # Read in bed file to dataframe
        binChr = np.loadtxt(self.binData, usecols = (0,), dtype = str)
        binStart = np.loadtxt(self.binData, usecols = (1,), dtype = np.uint32)
        binEnd = np.loadtxt(self.binData, usecols = (2,), dtype = np.uint32)
        binIndex = np.arange(0, len(binStart), 1, dtype = np.uint32)
        if len(binChr) != len(binStart) or len(binStart) != len(binEnd):
            raise IOError('Invalid bed file supplied')
        # Create output data variables
        binDict = collections.OrderedDict()
        binNames = []
        count = len(binStart)
        # Create arrays for every chromosome
        for chrom in np.unique(binChr):
            # Extract arrays for chromosome
            indices = binChr == chrom
            chrChr = binChr[indices]
            chrStart = binStart[indices]
            chrEnd = binEnd[indices]
            chrIndex = binIndex[indices]
            # Check for overlapping bins
            binArray = np.empty(2 * sum(indices), dtype = np.uint32)
            binArray[0::2] = chrStart
            binArray[1::2] = chrEnd
            for i in xrange(len(binArray) - 1):
                if binArray[i + 1] <= binArray[i]:
                    raise IOError('Bed file contains overlapping bins')
            # Store output data
            for bin in np.transpose([chrChr, chrStart, chrEnd]):
                binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
            binDict[chrom] = {'start' : chrStart, 'end' : chrEnd,
                'index' : chrIndex, 'count' : len(chrStart)}
        # Return data
        return(binDict, binNames, count)
    
    def findBinIndex(self, chrom, position):
        # Convert position format
        position = np.uint32(position)
        # Extract bin location
        try:
            binLocation =  self.binDict[chrom]['end'].searchsorted(position)
        except KeyError:
            return('nochr')
        # Extract bin data
        try:
            binStart = self.binDict[chrom]['start'][binLocation]
            binIndex = self.binDict[chrom]['index'][binLocation]
        except IndexError:
            return('nobin')
        # Extract bin name/number
        if binStart <= position:
            return(binIndex)
        else:
            return('nobin')
    
    def writeBed(self, fileName):
        with open(fileName, 'w') as outFile:
            for name in self.binNames:
                chrom, interval = name.split(':')
                start, end = interval.split('-')
                outFile.write('%s\t%s\t%s\n' %(chrom, start, end))
    
    def matrixProcess(self, inputQueue, outPipe):
        # Create matrix and log array
        matrix = np.zeros((self.binCount, self.binCount),
            dtype = np.uint32)
        logData = np.zeros(4, dtype = np.uint32)
        # Extract data from pipe 
        for fragPair in iter(inputQueue.get, None):
            # Count total
            logData[0] += 1
            # Strip and process line
            fragPair = fragPair.strip().split('\t')
            # Extract location of fragends
            indices = []
            fragends = [fragPair[0:2], fragPair[3:5]]
            for fragChr, fragLoc in fragends:
                index = self.findBinIndex(fragChr, fragLoc)
                if isinstance(index, np.uint32):
                    indices.append(index)
                else:
                    indices = index
                    break
            # Check that two bin  indexes have been identified
            if isinstance(indices, list):
                # Count accepted ligations and data to array
                logData[3] += 1
                matrix[indices[0]][indices[1]] += 1
                matrix[indices[1]][indices[0]] += 1
            # Count incorrect indices
            elif indices == 'nochr':
                logData[1] += 1
            elif indices == 'nobin':
                logData[2] += 1
            else:
                raise ValueError('unrecognised bin index')
        outPipe.send((matrix, logData))
    
    def generateMatrix(self, fragendFile, threads=1):
        # Manage thread number
        if threads > 2:
            threads -= 1
        # Open input file
        if fragendFile.endswith('.gz'):
            sp = subprocess.Popen(["zcat", fragendFile],
                stdout = subprocess.PIPE)
            fh = cStringIO.StringIO(sp.communicate()[0])
        else:
            fh = open(fragendFile)
        # Create queue
        fragQueue = multiprocessing.Queue()
        # Start processes to count interactions
        processData = []
        for _ in range(threads):
            # Create pipes and process
            pipeReceive, pipeSend = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = self.matrixProcess,
                args = (fragQueue, pipeSend)
            )
            process.start()
            pipeSend.close()
            # Strore process and pipe data
            processData.append((process,pipeReceive))
        # Add input data to queue
        for line in fh:
            fragQueue.put(line)
        # Add termination values to queue and close
        for _ in processData:
            fragQueue.put(None)
        fragQueue.close()
        # Extract data from processes and terminate
        for count, (process, pipe) in enumerate(processData):
            # Extract data and close process and pipes
            processMatrix, processLog = pipe.recv()
            process.join()
            pipe.close()
            # Add matrix count data
            if count:
                finalMatrix += processMatrix
                finalLog += processLog
            else:
                finalMatrix = processMatrix
                finalLog = processLog
        # Close input file
        if not fragendFile.endswith('.gz'):
            fh.close()
        # Return data
        return(finalMatrix, finalLog)

def normaliseMatrix(
        inMatrix, excludeBin, path = 'ic_mes', memory = 4000, maxIter = 1000
    ):
    ''' Function normalises HiC data using the iterative correction
    algorithm generated by Imakaev et al. The algorithm is incoporated
    into the Hi-Corrector package. Function takes 5 arguments:
    
    1)  inMatrix - either an string listing input matrix file or a numpy
        array of the raw counts. It is presumed that the input matrix file
        will have a header.
    2)  excludeBin - either an integer encoding the minimum bin count or a
        boolean numpy array with a value for each bin. A value of True
        will lead to a bin being excluded.
    3)  path - Path to ic_mes executable.
    4)  memory - memory size (MB)
    5)  maxIter - Maximum iterations
    
    '''
    # Read input file if supplied
    if isinstance(inMatrix, str):
        inMatrix = np.loadtxt(fname = inMatrix, dtype = np.uint32,
            delimiter = '\t', skiprows = 1)
    # Count bin no
    if inMatrix.shape[0] == inMatrix.shape[1]:
        binNo = inMatrix.shape[0]
    else:
        raise IOError('Matrix must be square')
    # Find bins to exclude from normalisation
    if isinstance(excludeBin, int):
        colSums = inMatrix.sum(axis=0)
        excludeBin = colSums < excludeBin
    if binNo != len(excludeBin):
        raise IOError('Exclude bin array must be equal to bin number')
    # Set excluded bins to zero
    inMatrix[excludeBin,:] = 0
    inMatrix[:,excludeBin] = 0
    # Create temporary files
    dirName = tempfile.mkdtemp()
    biasFile = os.path.join(dirName, 'bias.txt')
    tempMatrix = os.path.join(dirName, '.tempMatrix.txt')
    # Save altered matrix to file
    np.savetxt(tempMatrix, inMatrix, '%s', '\t')
    # Calculate bias 
    biasCommand = [path, tempMatrix, str(memory), str(binNo), str(maxIter),
        '0', '0', biasFile]
    biasLog = subprocess.check_output(biasCommand)
    # Extract bias val ues
    biasFactors = np.loadtxt(fname = biasFile, dtype = 'float',
        delimiter = '\t')
    # Replace zero values in biasFactors
    biasFactors[biasFactors == 0] = 1
    # Create bias matrix
    biasArray = biasFactors * biasFactors[:, np.newaxis]
    # Apply bias matrix
    normArray = inMatrix / biasArray
    # Make each column sum to 1
    normColSums = normArray.sum(axis = 0)
    normColSums[normColSums == 0] = 1
    normArray = normArray / normColSums
    # Remove temporary files
    os.remove(tempMatrix)
    os.remove(biasFile)
    os.removedirs(dirName)
    # Return data
    return(normArray, biasFactors)

def normaliseMatrixProcess(
        args
    ):
    # Extract parameters
    inMatrix, excludeBin, path, memory, maxIter = args
    # Extract bin names
    with gzip.open(inMatrix, 'r') as infile:
        binNames = infile.next().strip()
    # Generate output files
    outPrefix = inMatrix[:-15]
    nMatrix = outPrefix + '.normMatrix.gz'
    binData = outPrefix + '.binData'
    dstData = outPrefix + '.dist.gz'
    # Perform normalisation
    normMatrix, bias = normaliseMatrix(inMatrix = inMatrix,
        excludeBin = excludeBin, path = path, memory = memory,
        maxIter = maxIter)
    # Save normalised interactions
    np.savetxt(nMatrix, normMatrix, '%.6f', '\t',
        header = binNames, comments = '')
    # Extract bin level data
    maskMatrix = analyseInteraction.maskMatrix(nMatrix)
    maskMatrix.binDF['bias'] = bias
    maskMatrix.binDirection()
    maskMatrix.binDistance()
    maskMatrix.binDF.to_csv(binData, '\t', '')
    # Extract global data
    distance = maskMatrix.combinedDistance()
    np.savetxt(dstData, distance, '%s', '\t')
