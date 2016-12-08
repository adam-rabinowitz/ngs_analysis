import collections
import gzip
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
import itertools
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

class compare_paired_matrices(object):
    
    def __init__(
            self, samples1, samples2, conditions, indir,
            suffix = 'normMatrix.gz'
        ):
        # Check arguments
        if not isinstance(samples1, list) or len(samples1) == 0:
            raise TypeError('prefix1 must be a list of length >= 1')
        if not isinstance(samples2, list) or len(samples2) == 0:
            raise TypeError('prefix1 must be a list of length >= 1')
        if not os.path.isdir(indir):
            raise IOError('could not find indir')
        if (not isinstance(conditions, list)
            or len(conditions) != 2
            or not isinstance(conditions[0], str)
            or not isinstance(conditions[1], str)):
            raise TypeError('conditions must be a list of two strings')
        if not isinstance(suffix, str):
            raise TypeError('suffix must be a string')
        # Store supplied variables
        self.samples1 = samples1
        self.samples2 = samples2
        self.conditions = conditions
        self.indir = indir
        self.suffix = suffix
        # Create matrix dictionary and check pairings
        self.matrices = self.__create_matrix_dictionary()
        self.matrixnames = self.__check_matrix_pairing()
        
    def __create_matrix_dictionary(self):
        # Check that all prefixes are unique strings
        for s1, s2 in itertools.permutations(self.samples1 + self.samples2, 2):
            if not isinstance(s1, str) or not isinstance(s2, str):
                raise TypeError('prefixes must be lists of strings')
            if s1.startswith(s2):
                raise ValueError('{} & {} prefixes not unique'.format(s1, s2))
        # Create dictionary to store extracted matrices
        matrixDict = collections.OrderedDict()
        for condition, samples in zip(
                self.conditions, [self.samples1, self.samples2]
            ):
            matrixDict[condition] = collections.OrderedDict()
            for s in samples:
                matrixDict[condition][s] = []
        # Extract matrices for each prefix
        fileList = os.listdir(self.indir)
        fileList = [f for f in fileList if f.endswith(self.suffix)]
        for condition in matrixDict:
            for sample in matrixDict[condition]:
                for f in fileList:
                    if f.startswith(sample):
                        matrixDict[condition][sample].append(f)
        # Sort matrix lists and add full path
        for condition in matrixDict:
            for sample in matrixDict[condition]:
                matrixList = matrixDict[condition][sample]
                matrixList.sort()
                matrixList = [os.path.join(self.indir, m) for m in matrixList]
                matrixDict[condition][sample] = matrixList
        # Return data
        return(matrixDict)
    
    def __check_matrix_pairing(self):
        # Check that all matrices are paired for all samples
        reference = None
        for condition in self.matrices:
            for sample in self.matrices[condition]:
                # Extract matrix list and check length
                matrixList = self.matrices[condition][sample]
                if len(matrixList) == 0:
                    raise ValueError('No matrix found for {}'.format(sample))
                # Extract matrix names
                regx = re.compile('^.*?/{}([^/]+){}$'.format(sample, self.suffix))
                matrixNames = [regx.sub('\\1', x) for x in matrixList]
                matrixNames = [x.strip('.') for x in matrixNames]
                matrixNames.sort()
                # Check names are consitent
                if reference is None:
                    reference = matrixNames
                if matrixNames != reference:
                    print('Reference: {}'.format(', '.join(reference)))
                    print('Comparison: {}'.format(', '.join(matrixList)))
                    raise ValueError('Matrix names do not match')
        # Return data
        return(reference)

    def __prob_matrix(self, path):
        # Read in matrix
        matrix = np.loadtxt(path, delimiter='\t', skiprows=1,
            dtype=np.float64)
        # Check matrix is square
        m, n = matrix.shape
        if m != n:
            raise ValueError('{} is not square'.format(path))
        # Check matrix is symetrical to six decimal places
        if not np.allclose(matrix,  matrix.T, atol=1.01e-6, rtol=0):
            self.__extract_nonsymetrical_pairs(matrix)
        # Return matrix
        return(matrix)
    
    def __dist_matrix(self, path):
        # Extract bin names
        with gzip.open(path) as inFile:
            binNames = inFile.next().strip().split('\t')
        # Extract and check bin data
        binData = np.array([re.split('[:-]', x) for x in binNames])
        if not np.alltrue(binData.T[0] == binData.T[0][1]):
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

    def __extract_nonsymetrical_pairs(self, matrix, maxno=10):
        different = np.where(matrix != matrix.T)
        count = 0
        print('Non-symetrical values found. Examples follow:')
        for d1, d2 in zip(different[0], different[1]):
            output = '{}\t{}\t{}\t{}'.format(d1, d2, matrix[d1, d2], matrix[d2, d1])
            print(output)
            count += 1
            if count == maxno:
                break
        raise ValueError('Non symetrical matrix found')
    
    def extract_dist_prob(self):
        # Create output dataframe
        outDF = pd.DataFrame(
            columns = ['cond', 'repl', 'smpl', 'mtrx', 'dist', 'prob'])
        # Extract probabilities for input matrices
        for condition in self.matrices:
            for replicate, sample in enumerate(self.matrices[condition]):
                for path, name in zip(
                        self.matrices[condition][sample], self.matrixnames
                    ):
                    # Create matrices
                    probMatrix, distMatrix = self.__prob_dist_matrix(path)
                    # Extract data from lower triangles
                    trilIndices = np.tril_indices(probMatrix.shape[0])
                    probData = probMatrix[trilIndices]
                    distData = distMatrix[trilIndices]
                    # Create dataframe and concat to output
                    pathDF = pd.DataFrame({
                        'cond' : pd.Series([condition] * len(probData)),
                        'repl' : pd.Series([replicate + 1] * len(probData)),
                        'smpl' : pd.Series([sample] * len(probData)),
                        'mtrx' : pd.Series([name] * len(probData)),
                        'dist' : distData,
                        'prob' : probData
                    }, columns = ['cond', 'repl', 'smpl', 'mtrx', 'dist',
                        'prob'])
                    outDF = pd.concat((outDF, pathDF), axis=0)
        return(outDF)
    
    def mean_matrix_dist_prob(self, rmzero=True):
        # Check arguments:
        if not isinstance(rmzero, bool):
            raise TypeError('rmzero must be bool')
        # Extract probabilities and remove zeros, if requested
        probData = self.extract_dist_prob()
        if rmzero:
            probData = probData[probData['prob'] > 0]
        # Split the data and create output dataframe
        g = probData.groupby(['smpl', 'mtrx', 'dist'])
        outDF = pd.DataFrame(index = g.groups.keys(), columns=[
            'cond', 'repl', 'smpl', 'mtrx', 'dist', 'no', 'prob'])
        outDF = outDF.sort_index()
        # Populate dataframe
        for key, data in g:
            # Check all conditions and replicate data is identical
            if (data['repl'] != data['repl'].iloc[0]).any():
                raise ValueError('replicate not consistent across samples')
            if (data['cond'] != data['cond'].iloc[0]).any():
                raise ValueError('condition not consistent across samples')
            # Create and store output
            output = [data['cond'].iloc[0], data['repl'].iloc[0], key[0],
                key[1], key[2], data['prob'].size, data['prob'].mean()]
            outDF.loc[key] = output
        # Reindex and return dataframe
        outDF.index = np.arange(len(outDF))
        return(outDF)

    def calculate_dist_pvalue(self, rmzero=True, minvalues=10):
        # Extract distances for input matrices
        distProb = self.extract_dist_prob()
        splitDist = distProb.groupby('dist')
        # Create output columns
        colNames = []
        for condition in self.matrices:
            for sample in self.matrices[condition]:
                colNames.append('{}_{}_no'.format(condition, sample))
                colNames.append('{}_{}_mean'.format(condition, sample))
        for condition in self.matrices:
            colNames.append('{}_no'.format(condition))
            colNames.append('{}_mean'.format(condition))
        colNames.extend(['pvalue', 'fdr'])
        # Create output dataframe
        outDF = pd.DataFrame(
            columns = colNames, index = splitDist.groups.keys())
        outDF = outDF.sort_index()
        # Loop through data and calculate results
        for dist, data in splitDist:
            # Remove zero values
            if rmzero:
                data = data[data['prob'] > 0]
            # Extract data for conditions and samples
            condValues = []
            for cond in self.matrices:
                # Extract data for condition
                condData = data[data['cond'] == cond]
                condProb = condData['prob']
                condValues.append(condProb)
                # Add condition data to output
                colPrefix = '{}_'.format(cond)
                outDF.loc[dist, colPrefix + 'no'] = condProb.size
                outDF.loc[dist, colPrefix + 'mean'] = condProb.mean()
                for smpl in self.matrices[cond]:
                    # Extract data for sample
                    smplData = condData[condData['smpl'] == smpl]
                    smplProb = smplData['prob']
                    # Add sample data to output
                    colPrefix = '{}_{}_'.format(cond, smpl)
                    outDF.loc[dist, colPrefix + 'no'] = smplProb.size
                    outDF.loc[dist, colPrefix + 'mean'] = smplProb.mean()
            # Calculate pvalues
            prob1, prob2 = condValues
            if prob1.size >= minvalues and prob2.size >= minvalues:
                outDF.loc[dist, 'pvalue'] = mannwhitneyu(prob1, prob2)[1]
        # Sort data, add fdr and return
        pvalueIndex = outDF.index[~outDF['pvalue'].isnull()]
        outDF.loc[pvalueIndex, 'fdr'] = multipletests(
            outDF.loc[pvalueIndex, 'pvalue'], method='fdr_bh')[1]
        return(outDF)
    
    def extract_dist_quantile(
            self, quantile
        ):
        # Create output dataframe
        outDF = pd.DataFrame(columns=['cond', 'smpl', 'quan', 'dist'])
        # Loop through input queue
        for cond in self.matrices:
            for smpl in self.matrices[cond]:
                for path in self.matrices[cond][smpl]:
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
                        quantDF['cond'] = [cond] * len(quantile)
                        quantDF['smpl'] = [smpl] * len(quantile)
                        quantDF['quan'] = quantile
                        # Calculate quantile distances and store
                        for q in quantile:
                            subsetDF = distDF[distDF['cumsum'] >= q]
                            quantDF.loc[q, 'dist'] = subsetDF['dist'].iloc[0]
                        dfList.append(quantDF)
                    # Append results to output
                    outDF = pd.concat([outDF] + dfList, axis=0)
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
        # Create colnames for output dataframe
        colNames = []
        for condition in self.matrices:
            for sample in self.matrices[condition]:
                colNames.append('{}_{}_no'.format(condition, sample))
                colNames.append('{}_{}_mean'.format(condition, sample))
        for condition in self.matrices:
            colNames.append('{}_no'.format(condition))
            colNames.append('{}_mean'.format(condition))
        colNames.extend(['pvalue', 'fdr'])
        # Create output dataframe
        outDF = pd.DataFrame(index=quantile, columns=colNames)
        outDF = outDF.sort_index()
        # Extract quantile distance data
        quantData = self.extract_dist_quantile(quantile)
        splitQuant = quantData.groupby('quan')
        for q, data in splitQuant:
            # Extract data for conditions and samples
            condValues = []
            for cond in self.matrices:
                # Extract data for condition
                condData = data[data['cond'] == cond]
                condDist = condData['dist']
                condValues.append(condDist)
                # Add condition data to output
                colPrefix = '{}_'.format(cond)
                outDF.loc[q, colPrefix + 'no'] = condDist.size
                outDF.loc[q, colPrefix + 'mean'] = condDist.mean()
                for smpl in self.matrices[cond]:
                    # Extract data for sample
                    smplData = condData[condData['smpl'] == smpl]
                    smplDist = smplData['dist']
                    # Add sample data to output
                    colPrefix = '{}_{}_'.format(cond, smpl)
                    outDF.loc[q, colPrefix + 'no'] = smplDist.size
                    outDF.loc[q, colPrefix + 'mean'] = smplDist.mean()
            # Calculate pvalues
            dist1, dist2 = condValues
            if dist1.size >= minvalues and dist2.size >= minvalues:
                outDF.loc[q, 'pvalue'] = mannwhitneyu(dist1, dist2)[1]
        # Add fdr and return
        pvalueIndex = outDF.index[~outDF['pvalue'].isnull()]
        outDF.loc[pvalueIndex, 'fdr'] = multipletests(
            outDF.loc[pvalueIndex, 'pvalue'], method='fdr_bh')[1]
        return(outDF)
