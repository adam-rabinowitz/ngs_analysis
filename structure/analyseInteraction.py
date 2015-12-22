import re
import collections
import numpy as np
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
        # Extract chromsomal matrices and dataframes
        self.regionData = collections.OrderedDict()
        # Process regions if region data is provided
        if regions:
            print 'Regions not implemented'
        # Else seperate regions by chromosome
        else:
            for c in np.unique(self.binDF['chr']):
                # Extract data frame
                indices = np.where(self.binDF['chr'] == c)[0]
                regionDF = binDF.loc[indices]
                # Extract sub region
                probMatrix = self.matrix[indices,:][:,indices]
                # Create distance matrix
                middle = np.mean([regionDF['start'],regionDF['end']],index=0)
                distMatrix = np.abs(middle - middle[:, None])
                self.regionData[c] = (regionDF, probMatrix, distMatrix)
    
    def directionality(self):
        for r in self.regionData:
            
    
    def meanDirectionality(self):
        
        
        
