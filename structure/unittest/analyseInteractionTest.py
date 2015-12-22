import os
import unittest
import tempfile
import numpy as np
import pandas as pd
from ngs_analysis.structure import analyseInteraction


class InteractionTestCase(unittest.TestCase):
    
    def setUp(self):
        ''' Create data for unittest.'''
        # Create file names
        self.dirName = tempfile.mkdtemp()
        self.inMatrix = self.dirName + '/test.matrix'
        # Create matrix
        regions = ['chr1:1-10','chr1:11-41','chr1:42-49','chr2:1-21','chr2:22-30',
            'chr2:31-46','chr2:47-51','chr3:1-24']
        matrix = np.reshape(np.arange(len(regions)**2,dtype=float),
            (len(regions),len(regions)))
        np.savetxt(self.inMatrix, matrix, '%s', '\t', header = '\t'.join(regions),
            comments = '')
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        for file in [self.inMatrix]:
            if os.path.isfile(file):
                os.remove(file)
        os.removedirs(self.dirName)

class TestMatrixDivision(InteractionTestCase):

    def test_matrix_subdivision(self):
        ''' Check slicing of matrices and data frames '''
        subMatrix = analyseInteraction.subMatrix(self.inMatrix)
        # Check dataframe for chromosome 1
        df = pd.DataFrame()
        df['chr'] = np.array(['chr1']*3)
        df['start'] = np.array([1,11,42])
        df['end'] = np.array([10,41,49])
        self.assertTrue(all(df == subMatrix.regionData['chr1'][0]))
        # Check submatrix for chromosome 1
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr1'][1],
            np.array(
                [[0.,1.,2.],
                [8.,9.,10.],
                [16.,17.,18.]]
            )))
        # Check ditance matrix for chromosome 1
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr1'][2],
            np.array(
                [[0.,20.5,40.],
                [20.5,0.,19.5],
                [40.,19.5,0.]]
            )))
        # Check dataframe for chromosome 2
        df = pd.DataFrame()
        df['chr'] = np.array(['chr2']*4)
        df['start'] = np.array([1,22,31,47])
        df['end'] = np.array([21,30,46,51])
        df.index = [3,4,5,6]
        self.assertTrue(all(df == subMatrix.regionData['chr2'][0]))
        # Check submatrix for chromsome 2
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr2'][1],
            np.array(
                [[27.,28.,29.,30.],
                [35.,36.,37.,38.],
                [43.,44.,45.,46.],
                [51.,52.,53.,54.]]
            )))
        # Check ditance matrix for chromosome 2
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr2'][2],
            np.array(
                [[0.,15.,27.5,38.],
                [15.,0.,12.5,23.],
                [27.5,12.5,0.,10.5],
                [38.,23.,10.5,0.]]
            )))
        # Check dataframe for chromsome 3
        df = pd.DataFrame()
        df['chr'] = np.array(['chr3'])
        df['start'] = np.array([1])
        df['end'] = np.array([24])
        df.index = [7]
        self.assertTrue(all(df == subMatrix.regionData['chr3'][0]))
        # Check submatrix for chromosome 3
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr3'][1],
            np.array([
                [63.]
            ])))
        # Check ditance matrix for chromosome 3
        self.assertTrue(np.array_equal(
            subMatrix.regionData['chr3'][2],
            np.array([
                [0.]
            ])))


class TestMatrixDivision2(InteractionTestCase):

    def test_matrix_subdivision(self):
        ''' Check slicing of matrices and data frames '''
        subMatrix = analyseInteraction.subMatrix2(self.inMatrix)

suite = unittest.TestLoader().loadTestsFromTestCase(TestMatrixDivision)
unittest.TextTestRunner(verbosity=3).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(TestMatrixDivision2)
unittest.TextTestRunner(verbosity=3).run(suite)
