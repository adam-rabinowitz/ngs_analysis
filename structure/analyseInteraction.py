import re
import collections
import numpy as np
import numpy.ma as ma
import pandas as pd


class maskMatrix(object):
    
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

    def directionality(self):
        # Create series to store data
        self.binDF['self'] = self.binDF['inter'] = self.binDF['up'] = self.binDF['down'] = self.binDF['log2'] = np.nan
        # Loop through rows of the matrix
        for rowNo, row in enumerate(self.probMatrix):
            # Set none values if bin is entirely masked
            if ma.count(row) == 0:
                continue
            # Else calculate values
            else:
                # Extract self frequency
                selflig = row[rowNo]
                # Extract up frequency
                up = row[:rowNo]
                if ma.count(up) == 0:
                    up = 0.
                else:
                    up = up.sum()
                # Extract down frequency
                down = row[rowNo + 1:]
                if ma.count(down) == 0:
                    down = 0.
                else:
                    down = down.sum()
                # Calculate inter value
                inter = 1 - selflig - up - down
                # Calculate log2 value
                if up == 0:
                    if down == 0:
                        log2 = np.nan
                    else:
                        log2 = -np.inf
                elif down == 0:
                    log2 = np.inf
                else:
                    log2 = np.log2(up/down)
                # Store results
                self.binDF.loc[rowNo,['self','inter','up','down','log2']] = (
                    selflig, inter, up, down, log2)
