import collections
import gzip
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import mannwhitneyu

class analyse_interaction(object):
    
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

    def __distance_prob_generator(self, matrix):
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
    
    def __calc_quantile_metrics(
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
            for binName, binDF in self.__distance_prob_generator(matrix):
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
                target = self.__calc_quantile_metrics,
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

class compare_interaction(object):

    def __init__(self, matrixDict):
        # Store input dictionary and conditions
        if len(matrixDict) != 2:
            raise ValueError('Dictionary of two conditions must be supplied')
        # Check each condition has at least one matrix
        for paths in matrixDict.values():
            if not isinstance(paths, list):
                raise TypeError('Dictiuonary values must be lists')
            if len(paths) == 0:
                raise ValueError('Must be one or matrix for each condition')
            for p in paths:
                if not os.path.isfile(p):
                    raise IOError('File {} not found'.format(p))
        # Store data
        self.conditions = matrixDict.keys()
        self.matrixDict = matrixDict
    
    def __prob_matrix(self, path):
        # Read in matrix
        matrix = np.loadtxt(path, delimiter='\t', skiprows=1,
            dtype=np.float64)
        # Check matrix
        m, n = matrix.shape
        if m != n:
            raise ValueError('{} is not square'.format(path))
        if not np.alltrue(matrix == matrix.T):
            raise ValueError('{} is not symetrical'.format(path))
        # Return matrix
        return(matrix)
    
    def __dist_matrix(self, path):
        # Extract bin names
        with gzip.open(path) as inFile:
            binNames = inFile.next().strip().split('\t')
        # Extract and check bin data
        binData = np.array([re.split('[:-]', x) for x in binNames])
        if not np.alltrue(binData.T[0] == binData.T[0][1]):
            print(binData.T)
            raise ValueError('{} is not from single chromosome'.format(path))
        # Generate distance matrix
        centres = ((binData.T[1].astype(np.int64) +
            binData.T[2].astype(np.int64)) / 2).reshape(len(binData), 1)
        dist = abs(centres - centres.T)
        return(dist)
    
    def __prob_dist_matrix(self, path):
        probMatrix = self.__prob_matrix(path)
        distMatrix = self.__dist_matrix(path)
        if probMatrix.shape != distMatrix.shape:
            raise ValueError('{} has ambiguos bin numbers'.format('path'))
        return((probMatrix, distMatrix))
    
    def __extract_dist_prob(self):
        # Create output dataframe
        outDF = pd.DataFrame(columns = ['cond', 'dist', 'prob'])
        # Extract probabilities for input matrices
        for condition in self.matrixDict:
            for path in self.matrixDict[condition]:
                # Create matrices
                probMatrix, distMatrix = self.__prob_dist_matrix(path)
                # Extract data from lower triangles
                trilIndices = np.tril_indices(probMatrix.shape[0])
                probData = probMatrix[trilIndices]
                distData = distMatrix[trilIndices]
                # Create dataframe and concat to output
                pathDF = pd.DataFrame({
                    'cond' : pd.Series([condition] * len(probData)),
                    'dist' : distData,
                    'prob' : probData
                }, columns = ['cond', 'dist', 'prob'])
                outDF = pd.concat((outDF, pathDF), axis=0)
        return(outDF)
    
    def calculate_distance_difference(self, rmzero=True, minvalues=10):
        # Extract distances for input matrices
        distProb = self.__extract_dist_prob()
        splitDist = distProb.groupby('dist')
        # Create output dataframe
        colNames = [
            '{}_no'.format(self.conditions[0]),
            '{}_no'.format(self.conditions[1]),
            '{}_mean'.format(self.conditions[0]),
            '{}_mean'.format(self.conditions[1]),
            'pvalue',
            'fdr'
        ]
        outDF = pd.DataFrame(
            columns = colNames, index = splitDist.groups.keys())
        outDF = outDF.sort_index()
        # Loop through data and calculate results
        for dist, data in splitDist:
            # Remove zero values
            if rmzero:
                data = data[data['prob'] > 0]
            # Extract group data for first sample
            cond1 = self.conditions[0]
            prob1 = data['prob'][data['cond'] == cond1]
            outDF.loc[dist, cond1 + '_no'] = prob1.size
            outDF.loc[dist, cond1 + '_mean'] = prob1.mean()
            # Extract group data for second sample
            cond2 = self.conditions[1]
            prob2 = data['prob'][data['cond'] == cond2]
            outDF.loc[dist, cond2 + '_no'] = prob2.size
            outDF.loc[dist, cond2 + '_mean'] = prob2.mean()
            # Calculate p-value
            if prob1.size >= minvalues and prob2.size >= minvalues:
                outDF.loc[dist, 'pvalue'] = mannwhitneyu(prob1, prob2)[1]
        # Sort data, add fdr and return
        pvalueIndex = outDF.index[~outDF['pvalue'].isnull()]
        outDF.loc[pvalueIndex, 'fdr'] = multipletests(
            outDF.loc[pvalueIndex, 'pvalue'], method='fdr_bh')[1]
        return(outDF)
    
    def __extract_quantile_distance(
            self, quantile
        ):
        # Create output dataframe
        outDF = pd.DataFrame(columns=['cond', 'quan', 'dist'])
        # Loop through input queue
        for condition in self.matrixDict:
            for path in self.matrixDict[condition]:
                # Create matrices
                probMatrix, distMatrix = self.__prob_dist_matrix(path)
                dfList = []
                # Loop through columns
                for dist, prob in zip(distMatrix.T, probMatrix.T):
                    # Create dataframe listing distances
                    distDF = pd.DataFrame()
                    distDF['dist'] = dist
                    distDF['prob'] = prob
                    if distDF['prob'].sum() == 0:
                        continue
                    if not 1.05 > distDF['prob'].sum() > 0.95:
                        raise ValueError('Columns must add to 0 or ~1')
                    # Sort by distance and find cumulative frequence
                    distDF.sort_values('dist', inplace=True)
                    distDF['cumsum'] = distDF['prob'].cumsum()
                    # Create dataframe to store data
                    quantDF = pd.DataFrame(index = quantile)
                    quantDF['cond'] = [condition] * len(quantile)
                    quantDF['quan'] = quantile
                    # Calculate quantile distances and store
                    for q in quantile:
                        subsetDF = distDF[distDF['cumsum'] >= q]
                        quantDF.loc[q, 'dist'] = subsetDF['dist'].iloc[0]
                    dfList.append(quantDF)
                # Append results to output
                matrixDF = pd.concat(dfList, axis=0)
                outDF = pd.concat((outDF, matrixDF), axis=0)
        # Return data
        return(outDF)
    
    def calculate_quantile_pvalue(
            self, quantile, minvalues=10
        ):
        # Check arguments
        if isinstance(quantile, float):
            quantile = [quantile]
        elif isinstance(quantile, list):
            for q in quantile:
                if not isinstance(q, float):
                    raise TypeError('quantile list must contain floats')
        else:
            raise TypeError('quantile must be float or list of floats')
        # Create output dataframe
        colNames = [
            '{}_no'.format(self.conditions[0]),
            '{}_no'.format(self.conditions[1]),
            '{}_mean'.format(self.conditions[0]),
            '{}_mean'.format(self.conditions[1]),
            'pvalue',
            'fdr'
        ]
        outDF = pd.DataFrame(index=quantile, columns=colNames)
        outDF = outDF.sort_index()
        # Extract quantile distance data
        quantData = self.__extract_quantile_distance(quantile)
        splitQuant = quantData.groupby('quan')
        for q, data in splitQuant:
            # Extract group data for first sample
            cond1 = self.conditions[0]
            dist1 = data['dist'][data['cond'] == cond1]
            outDF.loc[q, cond1 + '_no'] = dist1.size
            outDF.loc[q, cond1 + '_mean'] = dist1.mean()
            # Extract group data for second sample
            cond2 = self.conditions[1]
            dist2 = data['dist'][data['cond'] == cond2]
            outDF.loc[q, cond2 + '_no'] = dist2.size
            outDF.loc[q, cond2 + '_mean'] = dist2.mean()
            # Calculate p-value
            if dist1.size >= minvalues and dist2.size >= minvalues:
                outDF.loc[q, 'pvalue'] = mannwhitneyu(dist1, dist2)[1]
        # Add fdr and return
        pvalueIndex = outDF.index[~outDF['pvalue'].isnull()]
        outDF.loc[pvalueIndex, 'fdr'] = multipletests(
            outDF.loc[pvalueIndex, 'pvalue'], method='fdr_bh')[1]
        return(outDF)



x = compare_interaction({
    'A' : [
        '/camp/stp/babs/working/rabinoa/yasu/mtrData/individualChrArmNormMatrices/NGS-13147.2000.III_ArmR.500.noself.normMatrix.gz',
        '/camp/stp/babs/working/rabinoa/yasu/mtrData/individualChrArmNormMatrices/NGS-13136.2000.III_ArmR.500.noself.normMatrix.gz'
    ],
    'B' : [
        '/camp/stp/babs/working/rabinoa/yasu/mtrData/individualChrArmNormMatrices/NGS-13146.2000.III_ArmR.500.noself.normMatrix.gz',
        '/camp/stp/babs/working/rabinoa/yasu/mtrData/individualChrArmNormMatrices/NGS-13138.2000.III_ArmR.500.noself.normMatrix.gz'
    ]
})
y = x.calculate_distance_difference()
print(y)










