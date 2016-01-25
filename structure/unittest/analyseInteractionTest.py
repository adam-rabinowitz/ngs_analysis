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
        self.inMatrix1 = self.dirName + '/test1.matrix'
        self.inMatrix2 = self.dirName + '/test2.matrix'
        # Create matrix 1
        regions1 = ['chr1:1-10','chr1:11-41','chr1:42-49','chr2:1-21',
            'chr2:22-30', 'chr2:31-46', 'chr2:47-51', 'chr2:52-63',
            'chr3:1-24']
        matrix1 = np.full((len(regions1),len(regions1)), 0.125)
        matrix1[4,:] = 0
        matrix1[:,4] = 0
        np.savetxt(self.inMatrix1, matrix1, '%s', '\t',
            header = '\t'.join(regions1), comments = '')
        # Create matrix 2
        regions2 = ['chr1:1-1', 'chr1:2-2', 'chr1:3-3', 'chr1:4-4', 'chr1:5-5']
        matrix2 = np.array([
            [0.2, 0.2, 0.1, 0.2, 0.3],
            [0.4, 0.1, 0.2, 0.2, 0.1],
            [0.1, 0.1, 0.4, 0.1, 0.3],
            [0.2, 0.2, 0.2, 0.2, 0.2],
            [0.1, 0.4, 0.1, 0.2, 0.2]
        ])
        np.savetxt(self.inMatrix2, matrix2, '%s', '\t',
            header = '\t'.join(regions2), comments = '')
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        for file in [self.inMatrix1, self.inMatrix2]:
            if os.path.isfile(file):
                os.remove(file)
        os.removedirs(self.dirName)

    def compna(self, s1, s2):
        ''' Create function to compare series containing nan values '''
        for v1, v2 in zip(s1, s2):
            if isinstance(v1, (str, int)):
                if v1 != v2:
                    return False
            elif np.isnan(v1):
                if not np.isnan(v2):
                    return False
            else:
                if v1 != v2:
                    return False
        return True

class TestInteractionAnalysis(InteractionTestCase):
    
    def test_dataframe_generation(self):
        ''' checking bin directionality calculation '''
        maskMatrix = analyseInteraction.maskMatrix(self.inMatrix1)
        maskMatrix.binDirection()
        df = maskMatrix.binDF
        self.assertTrue(self.compna(df['group'], pd.Series(
        ['chr1', 'chr1', 'chr1', 'chr2', np.nan, 'chr2', 'chr2', 'chr2', 'chr3'])))
        self.assertTrue(self.compna(df['up'], pd.Series(
            [0, 0.125, 0.25, 0, np.nan, 0.125, 0.25, 0.375, 0])))
        self.assertTrue(self.compna(df['down'], pd.Series(
            [0.25, 0.125, 0, 0.375, np.nan, 0.25, 0.125, 0, 0])))
        self.assertTrue(self.compna(df['inter'], pd.Series(
            [0.625, 0.625, 0.625, 0.5, np.nan, 0.5, 0.5, 0.5, 0.875])))
        self.assertTrue(self.compna(df['log2'], pd.Series(
            [np.nan, 0, np.nan, np.nan, np.nan, 0, 0, np.nan, np.nan])))

    def test_bin_distance1(self):
        ''' Test log2 calculation with standard window. '''
        maskMatrix = analyseInteraction.maskMatrix(self.inMatrix2)
        maskMatrix.binDirection()
        df = maskMatrix.binDF
        self.assertTrue(self.compna(df['log2'], pd.Series(
            [np.nan, 1, -1, 0])))
    
    def test_bin_distance2(self):
        ''' Test log2 calculation with short window. '''
        maskMatrix = analyseInteraction.maskMatrix(self.inMatrix2)
        maskMatrix.binDirection(maxl=1)
        df = maskMatrix.binDF
        self.assertTrue(self.compna(df['log2'], pd.Series(
            [np.nan, 1, 0, 0])))


suite = unittest.TestLoader().loadTestsFromTestCase(TestInteractionAnalysis)
unittest.TextTestRunner(verbosity=3).run(suite)

