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
        regions = ['chr1:1-10','chr1:11-41','chr2:1-21','chr2:22-30',
            'chr2:31-46','chr3:1-24']
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
        # Check output for chromosome 1
        self.assertTrue(np.array_equal(
            subMatrix.regMatrix['chr1'][0],
            np.array([
                [0.,1.],
                [6.,7.]
            ])))
        df = pd.DataFrame()
        df['chr'] = np.array(['chr1','chr1'])
        df['start'] = np.array([1,11])
        df['end'] = np.array([10,41])
        self.assertTrue(all(df == subMatrix.regMatrix['chr1'][1]))
        # Check output for chromosome 2
        self.assertTrue(np.array_equal(
            subMatrix.regMatrix['chr2'][0],
            np.array([
                [14.,15.,16.],
                [20.,21.,22.],
                [26.,27.,28.]
            ])))
        df = pd.DataFrame()
        df['chr'] = np.array(['chr2','chr2','chr2'])
        df['start'] = np.array([1,22,31])
        df['end'] = np.array([21,30,46])
        df.index = [2,3,4]
        self.assertTrue(all(df == subMatrix.regMatrix['chr2'][1]))
        # Check output for chromsome 3
        self.assertTrue(np.array_equal(
            subMatrix.regMatrix['chr3'][0],
            np.array([
                [35.]
            ])))
        df = pd.DataFrame()
        df['chr'] = np.array(['chr3'])
        df['start'] = np.array([1])
        df['end'] = np.array([24])
        df.index = [5]
        self.assertTrue(all(df == subMatrix.regMatrix['chr3'][1]))

suite = unittest.TestLoader().loadTestsFromTestCase(TestMatrixDivision)
unittest.TextTestRunner(verbosity=3).run(suite)
