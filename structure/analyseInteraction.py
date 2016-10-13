import collections
import gzip
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
from statsmodels.nonparametric.smoothers_lowess import lowess

class interaction_analysis(object):
    
    def __init__(self, matrixList):
        # Store matrix list and associated parameters
        self.matrixList = matrixList
        # Extract parameters and store regions
        regionDict = collections.defaultdict(list)
        for count, matrix in enumerate(matrixList):
            # Check existence of file
            if not os.path.isfile(matrix):
                raise IOError('Could not find input file')
            # Extract sample data
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
                    raise ValueError('Regions absent for some samples')
            else:
                self.regionList = regions

    def distance_prob_generator(self, matrix):
        # Extract bin names
        with gzip.open(matrix) as inFile:
            binNames = inFile.next().strip().split()
        # Create distance matrix
        start, end = zip(
            *[re.split(':|-', x)[1:] for x in binNames])
        centres = np.array([map(int, start), map(int, end)]).mean(axis=0)
        distMatrix = np.abs(centres - centres[:,None])
        # Read in matrix and remove columns
        probMatrix = np.loadtxt(
            matrix, dtype=np.float32, delimiter='\t', skiprows=1)
        for index, (prob, dist) in enumerate(zip(probMatrix.T, distMatrix.T)):
            binDF = pd.DataFrame()
            binDF['prob'] = prob
            binDF['dist'] = dist
            yield((binNames[index], binDF))
    
    def calc_quantile_metrics(
            self, inQueue, outQueue
        ):
        # Loop through input queue
        for matrix, quantile in iter(inQueue.get, None):
            # Extract bin data
            fileName = os.path.basename(matrix)
            sample, binSize, region, minCount = fileName.split('.')[:4]
            # Create output dataframe
            with gzip.open(matrix) as inFile:
                binNames = inFile.next().strip().split('\t')
            outDF = pd.DataFrame(index = binNames)
            outDF['sample'] = sample
            outDF['region'] = region
            outDF['bin'] = binNames
            # Create name for output series
            quantileName = quantile * 100
            if quantileName % 1:
                quantileName = 'Q' + str(quantileName)
            else:
                quantileName = 'Q' + str(int(quantileName))
            # Extract bin data
            for binName, binDF in self.distance_prob_generator(matrix):
                # Skip bins with low probabilites
                if binDF['prob'].sum() < quantile:
                    continue
                # Calculate median distance
                binDF['absdist'] = binDF['dist'].abs()
                binDF.sort_values('absdist', inplace=True)
                binDF['cumsum'] = binDF['prob'].cumsum()
                quantDist = binDF['dist'][binDF['cumsum'] >= quantile].iloc[0]
                # Add value to output df
                outDF.loc[binName, quantileName] = quantDist
            # Add output dataframe to out queue
            outQueue.put((sample, region, outDF))

    def calc_quantile(
            self, quantile, threads=1
        ):
        # Check arguments
        if not isinstance(quantile, float) or not 0 < quantile < 1:
            raise ValueError('quantile must be float between 0 and 1')
        # Create queues
        inQueue = multiprocessing.Queue()
        outQueue = multiprocessing.Queue()
        # Create processes
        processList = []
        for _ in range(threads):
            process = multiprocessing.Process(
                target = self.calc_quantile_metrics,
                args = (inQueue, outQueue)
            )
            process.start()
            processList.append(process)
        # Add data to queue
        for matrix in self.matrixList:
            inQueue.put((matrix, quantile))
        # Create ordered dictionary to store values
        output = collections.OrderedDict()
        for sample in self.sampleList:
            output[sample] = collections.OrderedDict()
            for region in self.regionList:
                output[sample][region] = None
        # Populate output dictionary with dataframes
        for _ in self.matrixList:
            sample, region, df = outQueue.get()
            output[sample][region] = df
        # Clean up
        for _ in range(threads):
            inQueue.put(None)
        for process in processList:
            process.join()
        # Concatenate regions and samples and return output
        for sample in output:
            output[sample] = pd.concat(
                output[sample].values(), axis = 0)
        output = pd.concat(
            output.values(), axis=0
        )
        output.index = np.arange(output.shape[0])
        return(output)

inDir = '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/'
inFiles = os.listdir(inDir)
matrixList = [os.path.join(inDir, x) for x in inFiles if x.endswith('.normMatrix.gz')]
x = interaction_analysis(matrixList)
y = x.calc_quantile(0.5, 6)
y.to_csv(os.path.join(inDir, 'Q50.df.txt'), sep='\t', index=False)

class maskMatrix(object):
    
    def __init__(self, matrix, regions='', overlap = False):
        # Extract bin names
        if matrix.endswith('.gz'):
            with gzip.open(matrix) as inFile:
                self.binNames = inFile.next().strip().split('\t')
        else:
            with open(matrix) as inFile:
                self.binNames = inFile.next().strip().split('\t')
        # Create bin dataframe
        binDF = pd.DataFrame()
        binDF['chr'], binDF['start'], binDF['end'] = zip(
            *[re.split(':|-', x) for x in self.binNames])
        binDF[['start', 'end']] = binDF[['start', 'end']].astype(int)
        binDF['chr'] = binDF['chr'].astype(str)
        binDF['name'] = self.binNames
        binDF['centre'] = np.mean([binDF['start'], binDF['end']], axis=0)
        self.binDF = binDF
        # Open probability matrix and check it is square
        probMatrix = np.loadtxt(matrix, skiprows = 1)
        if probMatrix.shape[0] != probMatrix.shape[1]:
            raise IOError('Matrix must be square')
        # Create mask matrix and remove low values
        chrArray = np.array(self.binDF['chr'])
        maskMatrix = chrArray != chrArray[:,None]
        lowValues = probMatrix.sum(axis=0) < 0.5
        maskMatrix[lowValues,:] = True
        maskMatrix[:,lowValues] = True
        # Add group data to dataframe
        groups = self.binDF['chr'].copy()
        groups[lowValues] = np.nan
        self.binDF['group'] = groups
        # Create masked probability and distance matrices
        self.probMatrix = ma.masked_array(probMatrix, mask = maskMatrix)
        centreArray = np.array(self.binDF['centre'])
        distMatrix = np.abs(centreArray - centreArray[:,None])
        self.distMatrix = ma.masked_array(distMatrix, mask = maskMatrix)
    
    def upDown(self, array, index):
        ''' Returns distance-paired upstream and downstream arrays.'''
        if index == 0:
            up = ma.masked_array([], mask=np.array([],dtype=bool))
        else:
            up = array[index-1::-1]
        down = array[index+1:]
        return(up, down)
    
    def unmaskedPair(self, a1, a2, maxl):
        ''' Returns indices of unmasked array pairs.'''
        maxl = min(len(a1), len(a2), maxl)
        if maxl == 0:
            return(np.array([]))
        masked = np.logical_or(a1.mask[:maxl], a2.mask[:maxl])
        indices = np.where(masked == False)[0]
        return(indices) 
    
    def binDirection(self, maxl = 10):
        ''' Extract interaction direction data for bins '''
        # Create dataframe to store data
        df = pd.DataFrame(columns = ['self', 'inter', 'up', 'down', 'log2'])
        self.binDF = pd.concat([self.binDF, df], axis=1)
        # Loop through rows of the matrix
        for rowNo, row in enumerate(self.probMatrix):
            # Set none values if bin is entirely masked
            if ma.count(row) == 0:
                continue
            # Else calculate values
            else:
                # Extract up and down arrays
                up, down = self.upDown(row, rowNo)
                # Extract probabilities
                selfProb = row[rowNo].sum()
                if up.count() == 0:
                    upProb = 0
                else:
                    upProb = up.sum()
                if down.count() == 0:
                    downProb = 0
                else:
                    downProb = down.sum()
                interProb = 1 - upProb - downProb - selfProb
                # Extract paired bins for log2 calculations
                indices = self.unmaskedPair(up, down, maxl)
                if len(indices) > 0:
                    # Calculate sum of paired up and down bins
                    upSum = up[indices].sum()
                    downSum = down[indices].sum()
                    # Calculate log2 ratio
                    if upSum == 0:
                        if downSum == 0:
                            log2 = np.nan
                        else:
                            log2 = -np.inf
                    elif downSum == 0:
                        log2 = np.inf
                    else:
                        log2 = np.log2(upSum/downSum)
                else:
                    log2 = np.nan
                # Store results
                self.binDF.loc[rowNo,['self','inter','up','down','log2']] = (
                    selfProb, interProb, upProb, downProb, log2)
    
    def binDistance(self):
        ''' Extract weighted mean interaction distance for bins '''
        # Create dataframe to store data
        df = pd.DataFrame(columns = ['dist'])
        self.binDF = pd.concat([self.binDF, df], axis=1)
        # Loop through rows of the matrix
        for rowNo, row in enumerate(self.distMatrix):
            # Set none values if bin is entirely masked
            if ma.count(row) == 0:
                continue
            # Else calculate values
            else:
                # Calculate distance and store results
                dist = ma.average(row, weights = self.probMatrix[rowNo])
                self.binDF.loc[rowNo,'dist'] = np.uint32(dist)
    
    def combinedDistance(self):
        ''' Extract lowess smooth interaction frequency for dataset '''
        # Extract probabilities
        prob = self.probMatrix[~self.probMatrix.mask]
        dist = self.distMatrix[~self.distMatrix.mask]
        # Return data
        return(np.array([dist, prob]).T)
