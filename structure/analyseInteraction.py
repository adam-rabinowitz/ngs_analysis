import re
import collections
import numpy as np
import numpy.ma as ma
import pandas as pd

class subMatrix(object):
    
    def __init__(self, matrix, regions='', overlap = False):
        # Load matrix and check shape
        self.matrix = np.loadtxt(matrix, skiprows = 1)
        width, height = self.matrix.shape
        if width != height:
            raise IOError('Matrix must be square')
        # Extract bin names
        if matrix.endswith('.gz'):
            with gzip.open(matrix) as inFile:
                self.binNames = inFile.readline().strip.split('\t')
        else:
            with open(matrix) as inFile:
                self.binNames = inFile.readline().strip().split('\t')
        # Create bin dataframe
        self.binDF = pd.DataFrame()
        self.binDF['chr'], self.binDF['start'], self.binDF['end'] = zip(
            *[re.split(':|-', x) for x in self.binNames])
        self.binDF[['start', 'end']] = self.binDF[['start', 'end']].astype(int)
        # Extract region data
        self.regionData = collections.OrderedDict()
        if regions:
            # Extract region data if regions supplied
            print 'Regions not implemented'
        else:
            # Extract region data by chromosome
            for c in np.unique(self.binDF['chr']):
                # Extract data frame
                indices = np.where(self.binDF['chr'] == c)[0]
                regionDF = self.binDF.loc[indices]
                # Extract sub region
                probMatrix = self.matrix[indices,:][:,indices]
                # Create distance matrix
                if len(regionDF.index) > 1:
                    middle = np.mean([regionDF['start'], regionDF['end']],
                        axis=0)
                    distMatrix = np.abs(middle - middle[:, None])
                else:
                    distMatrix = np.array([[0.]])
                self.regionData[c] = (regionDF, probMatrix, distMatrix)
    
    def directionality(self):
        for region in self.regionData:
            log2Data = np.empty(len(regionData[0].index))
            for rowNo, row in enumerate(region):
                if row.sum < 0.1:
                    log2Data[rowNo] = None
                else:
                    sumLeft = row[:rowNo].sum()
                    sumRight = row[rowNo + 1:].sum()
                    try:
                        ratio = sumLeft/sumRight
                    except ZeroDivisionError:
                        ratio = float('inf')
                    log2Data[rowNo] = np.log2(ratio)
            self.regionData[region][0]['log2Dir'] = log2Data

            
class subMatrix2(object):
    
    def __init__(self, matrix, regions='', overlap = False):
        # Extract bin names
        if matrix.endswith('.gz'):
            with gzip.open(matrix) as inFile:
                self.binNames = inFile.readline().strip.split('\t')
        else:
            with open(matrix) as inFile:
                self.binNames = inFile.readline().strip().split('\t')
        # Create bin dataframe
        binDF = pd.DataFrame()
        binDF['chr'], binDF['start'], binDF['end'] = zip(
            *[re.split(':|-', x) for x in self.binNames])
        binDF[['start', 'end']] = binDF[['start', 'end']].astype(int)
        binDF['chr'] = binDF['chr'].astype(str)
        binDF['centre'] =  np.mean([binDF['start'], binDF['end']], axis=0)
        self.binDF = binDF
        # Create mask matrix
        chrArray = np.array(self.binDF['chr'])
        maskMatrix = chrArray != chrArray[:,None]
        # Create masked probability matrix
        probMatrix = np.loadtxt(matrix, skiprows = 1)
        if probMatrix.shape[0] != probMatrix.shape[1]:
            raise IOError('Matrix must be square')
        self.probMatrix = ma.masked_array(probMatrix, mask = maskMatrix)
        # Create masked distance matrix
        centreArray = np.array(self.binDF['centre'])
        distMatrix = np.abs(centreArray - centreArray[:,None])
        self.distMatrix = ma.masked_array(probMatrix, mask = maskMatrix)
        print self.probMatrix
        print self.distMatrix
