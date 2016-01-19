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
            'chr2:31-46', 'chr2:47-51', 'chr2:52-63', 'chr3:1-24']
        matrix = np.full((len(regions),len(regions)), 0.125)
        matrix[5,:] = 0
        matrix[:,5] = 0
        np.savetxt(self.inMatrix, matrix, '%s', '\t', header = '\t'.join(regions),
            comments = '')
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        for file in [self.inMatrix]:
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
    
    def test_bin_directionality(self):
        ''' checking bin directionality calculation '''
        maskMatrix = analyseInteraction.maskMatrix(self.inMatrix)
        maskMatrix.directionality()
        df = maskMatrix.binDF
        self.assertTrue(self.compna(df['group'], pd.Series(
        ['chr1', 'chr1', 'chr1', 'chr2', 'chr2', np.nan, 'chr2', 'chr2', 'chr3'])))
        self.assertTrue(self.compna(df['up'], pd.Series(
            [0, 0.125, 0.25, 0, 0.125, np.nan, 0.25, 0.375, 0])))
        self.assertTrue(self.compna(df['down'], pd.Series(
            [0.25, 0.125, 0, 0.375, 0.25, np.nan, 0.125, 0, 0])))
        self.assertTrue(self.compna(df['inter'], pd.Series(
            [0.625, 0.625, 0.625, 0.5, 0.5, np.nan, 0.5, 0.5, 0.875])))
        self.assertTrue(self.compna(df['log2'], pd.Series(
            [-np.inf, 0, np.inf, -np.inf, -1, np.nan, 1, np.inf, np.nan])))

    def test_bin_distance(self):
        maskMatrix = analyseInteraction.maskMatrix(self.inMatrix)
        maskMatrix.directionality()
        maskMatrix.distance()
        maskMatrix.lowessDistance()


suite = unittest.TestLoader().loadTestsFromTestCase(TestInteractionAnalysis)
unittest.TextTestRunner(verbosity=3).run(suite)

