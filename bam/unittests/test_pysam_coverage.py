import numpy as np
import os
import pandas as pd
import unittest
from ngs_python.bam import pysam_coverage

class TestMeanCoverage(unittest.TestCase):
    
    def setUp(self):
        dirPath = os.path.dirname(os.path.realpath(__file__))
        bamFile = os.path.join(dirPath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamFile)
        self.intervals = [
            ('ref', 0, 2),
            ('ref', 5, 9),
            ('ref', 9, 13),
            ('ref', 19, 27),
            ('ref', 43, 49),
            ('ref', 48, 56),
            ('ref', 50, 54),
            ('ref', 64, 70)
        ]
        self.names = ['{}:{}-{}'.format(*x) for x in self.intervals]
    
    def test_all_reads(self):
        meanCov = self.cov.mean_coverage(self.intervals)
        exptCov = pd.Series([0, 1.25, 2, 2.5, 2.333, 0.5, 0, 0],
            index=self.names)
        self.assertTrue(np.allclose(meanCov, exptCov, rtol=0, atol=0.001))
    
    def test_duplicate_reads(self):
        meanCov = self.cov.mean_coverage(self.intervals, remove_dup=True)
        exptCov = pd.Series([0, 1.25, 2, 2.25, 2, 0.5, 0, 0],
            index=self.names)
        self.assertTrue(np.allclose(meanCov, exptCov, rtol=0, atol=0.001))
    
    def test_secondary_reads(self):
        meanCov = self.cov.mean_coverage(self.intervals, remove_secondary=True)
        exptCov = pd.Series([0, 0.5, 1, 1.625, 2.333, 0.5, 0, 0],
            index=self.names)
        self.assertTrue(np.allclose(meanCov, exptCov, rtol=0, atol=0.001))
    
    def test_read_quality_10(self):
        meanCov = self.cov.mean_coverage(self.intervals, map_quality=20)
        exptCov = pd.Series([0, 1.25, 2, 2.5, 1.333, 0, 0, 0],
            index=self.names)
        self.assertTrue(np.allclose(meanCov, exptCov, rtol=0, atol=0.001))
    
    def test_read_quality_20(self):
        meanCov = self.cov.mean_coverage(self.intervals, map_quality=10)
        exptCov = pd.Series([0, 1.25, 2, 2.5, 2.333, 0.25, 0, 0],
            index=self.names)
        self.assertTrue(np.allclose(meanCov, exptCov, rtol=0, atol=0.001))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMeanCoverage)
    unittest.TextTestRunner(verbosity=2).run(suite)
