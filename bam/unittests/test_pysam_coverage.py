import numpy as np
import os
import unittest
from ngs_python.bam import pysam_coverage

class test_mean_coverage_each(unittest.TestCase):
    
    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
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
    
    def test_all_reads(self):
        meancov = self.cov.mean_coverage_each(self.intervals)
        exptcov = np.array([0, 1.25, 2, 2.5, 2.333, 0.5, 0, 0])
        self.assertTrue(np.allclose(meancov, exptcov, rtol=0, atol=0.001))
    
    def test_duplicate_reads(self):
        meancov = self.cov.mean_coverage_each(self.intervals, remove_dup=True)
        exptcov = np.array([0, 1.25, 2, 2.25, 2, 0.5, 0, 0])
        self.assertTrue(np.allclose(meancov, exptcov, rtol=0, atol=0.001))
    
    def test_secondary_reads(self):
        meancov = self.cov.mean_coverage_each(self.intervals,
            remove_secondary=True)
        exptcov = np.array([0, 0.5, 1, 1.625, 2.333, 0.5, 0, 0])
        self.assertTrue(np.allclose(meancov, exptcov, rtol=0, atol=0.001))
    
    def test_read_quality_10(self):
        meancov = self.cov.mean_coverage_each(self.intervals, map_quality=20)
        exptcov = np.array([0, 1.25, 2, 2.5, 1.333, 0, 0, 0])
        self.assertTrue(np.allclose(meancov, exptcov, rtol=0, atol=0.001))
    
    def test_read_quality_20(self):
        meancov = self.cov.mean_coverage_each(self.intervals, map_quality=10)
        exptcov = np.array([0, 1.25, 2, 2.5, 2.333, 0.25, 0, 0])
        self.assertTrue(np.allclose(meancov, exptcov, rtol=0, atol=0.001))

class test_coverage_count_all(unittest.TestCase):
    
    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
    
    def test_all_reads(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)])
        expDict = {0:12, 1:20, 2:10, 3:17, 4:11}
        self.assertTrue(covDict == expDict)
    
    def test_duplicate_reads(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)],
            remove_dup=True)
        expDict = {0:12, 1:20, 2:19, 3:19}
        self.assertTrue(covDict == expDict)
    
    def test_secondary_reads(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)],
            remove_secondary=True)
        expDict = {0:15, 1:25, 2:10, 3:10, 4:10}
        self.assertTrue(covDict == expDict)
    
    def test_read_quality_10(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)],
            map_quality=10)
        expDict = {0:22, 1:10, 2:10, 3:17, 4:11}
        self.assertTrue(covDict == expDict)
    
    def test_read_quality_20(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)],
            map_quality=20)
        expDict = {0:25, 1:9, 2:13, 3:22, 4:1}
        self.assertTrue(covDict == expDict)
    
    def test_filter_all(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)],
            map_quality=20, remove_secondary=True, remove_dup=True)
        expDict = {0:28, 1:19, 2:23}
        self.assertTrue(covDict == expDict)
    
    def test_max_cov(self):
        covDict = self.cov.coverage_count_all([('ref', 0, 70)], max_cov=3)
        expDict = {0:12, 1:20, 2:10, 3:28}
        self.assertTrue(covDict == expDict)
    
    def test_multi_region(self):
        covDict = self.cov.coverage_count_all([('ref', 10, 20),
            ('ref', 30, 40), ('ref', 50, 60)])
        expDict = {0:4, 1:6, 2:2, 3:13, 4:5}
        self.assertTrue(covDict == expDict)

class test_coverage_histogram(unittest.TestCase):
    
    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
    
    def test_all_reads(self):
        covhist = self.cov.coverage_histogram([('ref', 0, 70)], max_cov=5)
        exphist = np.array([1, 0.828, 0.542, 0.4, 0.157, 0])
        self.assertTrue(np.allclose(covhist, exphist, rtol=0, atol=0.001))
    
    def test_filter_all(self):
        covhist = self.cov.coverage_histogram([('ref', 0, 70)], max_cov=5,
            map_quality=20, remove_secondary=True, remove_dup=True)
        exphist = np.array([1, 0.6, 0.328, 0, 0, 0])
        self.assertTrue(np.allclose(covhist, exphist, rtol=0, atol=0.001))
    
    def test_max_cov(self):
        covhist = self.cov.coverage_histogram([('ref', 0, 70)], max_cov=3)
        exphist = np.array([1, 0.828, 0.542, 0.4])
        self.assertTrue(np.allclose(covhist, exphist, rtol=0, atol=0.001))

class test_create_bins(unittest.TestCase):
    
    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
    
    def test_equal_multiple(self):
        bindict = self.cov.create_bins(10, True)
        expstart = np.array([0, 10, 20, 30, 40, 50, 60])
        expend = np.array([10, 20, 30, 40, 50, 60, 70])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))
        
    def test_equal_remainder(self):
        bindict = self.cov.create_bins(9, True)
        expstart = np.array([3, 12, 21, 30, 39, 48, 57])
        expend = np.array([12, 21, 30, 39, 48, 57, 66])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))
    
    def test_equal_toobig(self):
        bindict = self.cov.create_bins(71, True)
        self.assertEqual(bindict['ref']['start'].size, 0)
        self.assertEqual(bindict['ref']['end'].size, 0)
    
    def test_unequal_multiple(self):
        bindict = self.cov.create_bins(10, False)
        expstart = np.array([0, 10, 20, 30, 40, 50, 60])
        expend = np.array([10, 20, 30, 40, 50, 60, 70])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))
        
    def test_unequal_remainder1(self):
        bindict = self.cov.create_bins(9, False)
        expstart = np.array([0, 9, 18, 27, 36, 45, 54, 62])
        expend = np.array([9, 18, 27, 36, 45, 54, 62, 70])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))
        
    def test_unequal_remainder2(self):
        bindict = self.cov.create_bins(20, False)
        expstart = np.array([0, 18, 36, 53])
        expend = np.array([18, 36, 53, 70])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))
    
    def test_unequal_toobig(self):
        bindict = self.cov.create_bins(71, False)
        expstart = np.array([0])
        expend = np.array([70])
        self.assertTrue(np.all(bindict['ref']['start'] == expstart))
        self.assertTrue(np.all(bindict['ref']['end'] == expend))

class test_add_bin_count(unittest.TestCase):
    
    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
    
    def test_outside_5prime(self):
        bindict = self.cov.create_bins(9, True)
        self.cov.add_bin_count(bindict, 'ref', 2)
        expcount = np.array([0, 0, 0, 0, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
        
    def test_outside_3prime(self):
        bindict = self.cov.create_bins(9, True)
        self.cov.add_bin_count(bindict, 'ref', 66)
        expcount = np.array([0, 0, 0, 0, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))

    def test_boundary_5prime(self):
        bindict = self.cov.create_bins(9, True)
        self.cov.add_bin_count(bindict, 'ref', 29)
        expcount = np.array([0, 0, 1, 0, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
        
    def test_boundary_3prime(self):
        bindict = self.cov.create_bins(9, True)
        self.cov.add_bin_count(bindict, 'ref', 30)
        expcount = np.array([0, 0, 0, 1, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
       
class test_count_bin_overlaps(unittest.TestCase):

    def setUp(self):
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bamfile = os.path.join(dirpath, 'test_coverage.bam')
        self.cov = pysam_coverage.single_coverage(bamfile)
    
    def test_all_counts_refstart(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='refstart', rmDup=False,
            rmSec=False, rmSup=False)
        expcount = np.array([2, 1, 3, 1, 0, 1, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_readstart(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='readstart',
            rmDup=False, rmSec=False, rmSup=False)
        expcount = np.array([2, 1, 3, 1, 0, 0, 1])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_refend(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='refend',
            rmDup=False, rmSec=False, rmSup=False)
        expcount = np.array([0, 1, 2, 0, 4, 0, 1])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_readend(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='readend',
            rmDup=False, rmSec=False, rmSup=False)
        expcount = np.array([0, 1, 2, 0, 4, 1, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_mapq10(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=10, overlap='refstart',
            rmDup=False, rmSec=False, rmSup=False)
        expcount = np.array([2, 1, 3, 1, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_mapq20(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=20, overlap='refstart',
            rmDup=False, rmSec=False, rmSup=False)
        expcount = np.array([2, 1, 3, 0, 0, 0, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_rmdup(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='refstart', rmDup=True,
            rmSec=False, rmSup=False)
        expcount = np.array([2, 1, 2, 1, 0, 1, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_rmsec(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='refstart', rmDup=False,
            rmSec=True, rmSup=False)
        expcount = np.array([1, 1, 3, 1, 0, 1, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))
    
    def test_all_counts_rmsup(self):
        bindict = self.cov.count_bin_overlaps(
            binSize=10, binEqual=True, mapQ=0, overlap='refstart', rmDup=False,
            rmSec=False, rmSup=True)
        expcount = np.array([1, 1, 3, 1, 0, 1, 0])
        self.assertTrue(np.all(bindict['ref']['count'] == expcount))

if __name__ == '__main__':
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_mean_coverage_each)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_coverage_count_all)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_coverage_histogram)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_create_bins)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_add_bin_count)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(test_count_bin_overlaps)
    unittest.TextTestRunner(verbosity=2).run(suite)
