import collections
import gzip
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
import itertools

class tad_analysis(object):
    
    def __init__(self, matrixList):
        # Store matrix list and associated parameters
        self.matrixList = matrixList
        # Extract parameters and store regions
        regionDict = collections.defaultdict(list)
        for count, matrix in enumerate(matrixList):
            fileName = os.path.basename(matrix)
            if not fileName.endswith('.normMatrix.gz'):
                raise ValueError('Unexpected input files')
            sample, binSize, region, minCount = fileName.split('.')[:4]
            regionDict[sample].append(region)
            if count:
                if (not int(binSize) == self.binSize
                    or not int(minCount) == self.minCount):
                    raise ValueError('Sample parametes are different')
            else:
                self.binSize = int(binSize)
                self.minCount = int(minCount)
        # Check regions and store
        self.sampleList = regionDict.keys()
        self.sampleList.sort()
        for count, sample in enumerate(self.sampleList):
            regions = regionDict[sample]
            regions.sort()
            if count:
                if not regions == self.regionList:
                    print regionDict
                    print self.regionList
                    print regions
                    raise ValueError('Regions absent for some samples')
            else:
                self.regionList = regions


    def paired_prob_generator(self, matrix, max_length=100):
        # Extract bin names
        with gzip.open(matrix) as inFile:
            binNames = inFile.next().strip().split()
        # Read in matrix and remove columns
        prob = np.loadtxt(matrix, dtype=np.float32, delimiter='\t', skiprows=1)
        for index, column in enumerate(prob.T):
            # Extract upstream and downstream values
            upstream = column[:index]
            upstream = upstream[::-1]
            downstream = column[index + 1:]
            # Combined values into array
            length = min(upstream.shape[0], downstream.shape[0], max_length)
            pairedArray = np.vstack((upstream[:length], downstream[:length]))
            yield((binNames[index], pairedArray))
    
    def calc_log2_process(
            self, inQueue, outQueue
        ):
        # Loop through input queue
        for matrix, min_length, max_length in iter(inQueue.get, None):
            # Extract bin names
            with gzip.open(matrix) as inFile:
                binNames = inFile.next().strip().split()
            binArray = np.array([re.split('[:-]', x) for x in binNames])
            # Extract bin data
            fileName = os.path.basename(matrix)
            sample, binSize, region, minCount = fileName.split('.')[:4]
            # Create output dataframe
            outDF = pd.DataFrame(index = binNames)
            outDF['region'] = region
            outDF[sample] = np.nan
            # Loop through bins and extract pairs
            for binName, pairedArray in self.paired_prob_generator(
                    matrix, max_length
                ):
                pairedArray[pairedArray == 0] = np.nan
                ratio = pairedArray[0] / pairedArray[1]
                ratio = ratio[~np.isnan(ratio)]
                if ratio.shape[0] < min_length:
                    continue
                log2 = np.log2(ratio)
                mean = log2.mean()
                outDF.loc[binName, sample] = mean
            # Add output dataframe to out queue
            outQueue.put((sample, region, outDF))

    def calc_log2(
            self, min_length, max_length, threads
        ):
        # Check arguments
        if not isinstance(min_length, int) or min_length < 1:
            raise ValueError('max_length must be positive integer')
        if not isinstance(max_length, int) or max_length < 1:
            raise ValueError('max_length must be positive integer')
        if min_length > max_length:
            raise ValueError('min_length is higher than max_length')
        # Create queues
        inQueue = multiprocessing.Queue()
        outQueue = multiprocessing.Queue()
        # Create processes
        processList = []
        for _ in range(threads):
            process = multiprocessing.Process(
                target = self.calc_log2_process,
                args = (inQueue, outQueue)
            )
            process.start()
            processList.append(process)
        # Add data to queue
        for matrix in self.matrixList:
            inQueue.put((matrix, min_length, max_length))
        # Create ordered dictionary to store values
        outputDict = collections.OrderedDict()
        for sample in self.sampleList:
            outputDict[sample] = collections.OrderedDict()
            for region in self.regionList:
                outputDict[sample][region] = None
        # Populate output dictionary with dataframes
        for _ in self.matrixList:
            sample, region, df = outQueue.get()
            outputDict[sample][region] = df
        # Clean up
        for _ in range(threads):
            inQueue.put(None)
        for process in processList:
            process.join()
        # Concatenate regions for each sample
        for sample in outputDict:
            outputDict[sample] = pd.concat(
                outputDict[sample].values(), axis = 0)
            print(outputDict[sample].shape)
        # Create output dataframe
        binNames = outputDict[self.sampleList[0]].index
        binArray = np.array([re.split('[:|-]', x) for x in binNames])
        outputDF = pd.DataFrame(index = binNames)
        outputDF['bin'] = binNames
        outputDF['chr'] = binArray.T[0]
        outputDF['start'] = binArray.T[1]
        outputDF['end'] = binArray.T[2]
        outputDF['region'] = outputDict[self.sampleList[0]]['region']
        # Add sample data to dataframe and return
        for sample in outputDict:
            sampleDF = outputDict.pop(sample)
            outputDF[sample] = sampleDF[sample]
        return(outputDF)

#inputFiles = [
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8177.2000.III_ArmL.500.normMatrix.gz',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8178.2000.III_ArmL.500.normMatrix.gz',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8179.2000.III_ArmL.500.normMatrix.gz',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8177.2000.III_ArmR.500.normMatrix.gz',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8178.2000.III_ArmR.500.normMatrix.gz',
#    '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/NGS-8179.2000.III_ArmR.500.normMatrix.gz'
#]
inDir = '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/smc4G1Mitosis/'
inputFiles = os.listdir(inDir)
inputFiles = [os.path.join(inDir, x) for x in inputFiles if x.endswith('normMatrix.gz')]
ta = tad_analysis(inputFiles)
print(ta.regionList)
print(ta.sampleList)
#print(ta.binSize, ta.minCount)
df = ta.calc_log2(2, 100, 4)
print(df.shape)
df.to_csv('/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/smc4G1Mitosis/log2.df.txt', sep='\t', index=False)
