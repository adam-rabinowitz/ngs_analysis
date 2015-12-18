import collections
import multiprocessing
import os
import tempfile
import subprocess
import gzip
import numpy as np
import pandas as pd

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
                binEnd = np.arange(end, int(chrLength) + 1, maxWidth)
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
                binEnd = np.cumsum(binWidth)
                binStart = (binEnd - binWidth) + 1
            # Generate arrays containing bin index and chromosome name
            binIndex = np.arange(count, count + len(binStart))
            count += len(binStart)
            binChr = np.array([chrName] *  len(binStart))
            # Generate output data
            chrDF = pd.DataFrame({'chr': binChr, 'start':binStart,
                'end':binEnd, 'index':binIndex})
            chrDF = chrDF[['chr','start','end','index']]
            for bin in np.transpose([binChr, binStart, binEnd]):
                binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
            # Append output data
            binDict[chrName] = chrDF
        # Return dataframe
        return(binDict, binNames, count)

    def readBed(self):
        # Read in bed file to dataframe
        df = pd.read_csv(self.binData, sep = '\t', header = None)
        df = df.iloc[:,0:3]
        df.columns = ['chr','start','end']
        # Count bins and add index
        count = len(df.index)
        df['index'] = np.arange(0,count)
        # Extract bin names
        binNames = []
        for bin in np.transpose([df['chr'], df['start'], df['end']]):
            binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
        # Split data frames
        chrom = np.unique(df['chr'])
        binDF = collections.OrderedDict()
        for c in chrom:
            binDF[c] = df.loc[df['chr'] == c]
            binDF[c].index = np.arange(0, len(binDF[c].index))
        # Check that the bins are non-overlapping and ordered
        for chrDF in binDF.values():
            binArray = np.empty(2 * len(chrDF.index), dtype=int)
            binArray[0::2] = chrDF['start']
            binArray[1::2] = chrDF['end']
            for i in xrange(len(binArray) - 1):
                if binArray[i + 1] <= binArray[i]:
                    raise IOError('Bed file may contain overlapping bins')
        # Return data
        return(binDF, binNames, count)

    def findBinIndex(self, chrom, position):
        # Extract bin location
        try:
            binLocation =  np.searchsorted(
                self.binDict[chrom]['end'], position)
        except KeyError:
            return('nochr')
        # Extract bin data
        try:
            binStart, binIndex = self.binDict[chrom].loc[
                int(binLocation),['start','index']]
        except KeyError:
            return('nobin')
        # Extract bin name/number
        if binStart <= position:
            return(binIndex)
        else:
            return('nobin')
    
    def writeBed(self, fileName):
        outDF = pd.concat(self.binDict.values())
        np.savetxt(fileName, outDF, '%s', '\t')

def matrixProcess(genomeBin, inputQueue, outPipe):
    # Create matrix and log array
    matrix = np.zeros((genomeBin.binCount, genomeBin.binCount),
        dtype='int32')
    logData = np.zeros(4, dtype = 'int32')
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
            index = genomeBin.findBinIndex(fragChr, int(fragLoc))
            if isinstance(index, int):
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

def generateMatrix(fragendFile, genomeBin, threads=1):
    # Open input file
    if fragendFile.endswith('.gz'):
        inFile = gzip.open(fragendFile)
    else:
        inFile = open(fragendFile)
    # Create queue
    fragQueue = multiprocessing.Queue()
    # Start processes to count interactions
    processData = []
    for _ in range(threads):
        # Create pipes and process
        pipeReceive, pipeSend = multiprocessing.Pipe(False)
        process = multiprocessing.Process(
            target = matrixProcess,
            args = (genomeBin, fragQueue, pipeSend)
        )
        process.start()
        pipeSend.close()
        # Strore process and pipe data
        processData.append((process,pipeReceive))
    # Add input data to queue
    for line in inFile:
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
    # Return data
    return(finalMatrix, finalLog)

def normaliseMatrix(inMatrix, minCount):
    ''' Function normalises HiC data using the iterative correction algorithm
    generated by Imakaev et al. The algorithm is incoporated into the Hi-Corrector
    package. Function takes 3 arguments:
    1)  inFile - name of tab delimited input count matrix
    2)  outPrefix - prefix for the two output files '.bias' and '.normMatrix'
    3)  minCount - minimum count for a bin to be included in the analysis '''
    # Create temporary director
    dirName = tempfile.mkdtemp()
    # Create output file names
    biasFile = dirName + '/bias.txt'
    tempMatrix = dirName + '/tempMatrix.txt'
    # Count bin no
    binNo = inMatrix.shape[0]
    # Replace data from bins with less than minimum counts
    colSums = inMatrix.sum(axis=0)
    lowIndex = np.where(colSums < minCount)
    inMatrix[lowIndex] = 0
    inMatrix[:,lowIndex] = 0
    # Save processed counts to file
    np.savetxt(tempMatrix, inMatrix, '%s', '\t')
    # Calculate bias 
    biasCommand = ['/farm/babs/redhat6/software/Hi-Corrector1.1/bin/ic_mes', 
        tempMatrix, '5000', str(binNo), '1000', '0', '0', biasFile]
    biasLog = subprocess.check_output(biasCommand)
    # Extract bias values
    biasFactors = np.loadtxt(
        fname = biasFile,
        dtype = 'float',
        delimiter = '\t'
    )
    # Replace zero values in biasFactors
    biasFactors[biasFactors == 0] = 1
    # Create bias matrix
    biasArray = biasFactors[:, np.newaxis] * biasFactors
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
