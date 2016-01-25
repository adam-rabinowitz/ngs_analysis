import unittest
import collections
import tempfile
import os
import numpy as np
from ngs_python.structure import interactionMatrix

class InteractionTestCase(unittest.TestCase):
   
    def setUp(self):
        ''' Create data for unittest.'''
        # Create file names
        self.dirName = tempfile.mkdtemp()
        self.chrFile = self.dirName + '/test.chr'
        self.inBed1 = self.dirName + '/testIn1.bed'
        self.inBed2 = self.dirName + '/testIn2.bed'
        self.inFrag = self.dirName + '/testIn.frag'
        # Create files
        with open(self.chrFile, 'w') as chrFile:
            chrFile.write('%s\t%s\n%s\t%s' %('chr1',46,'chr2',33))
        with open(self.inBed1, 'w') as inBed:
            inBed.write('%s\t%s\t%s\n%s\t%s\t%s\n%s\t%s\t%s' %(
                'chr1',10,20,'chr1',30,40,'chr2',5,15))
        with open(self.inBed2, 'w') as inBed:
            inBed.write('%s\t%s\t%s\n%s\t%s\t%s\n%s\t%s\t%s' %(
                'chr1',10,20,'chr1',30,40,'chr1',40,41))
        with open(self.inFrag, 'w') as inFrag:
            inFrag.write('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' %(
                    'chr1\t4\t+\tchr1\t13\t-',
                    'chr2\t2\t+\tchr1\t33\t-',
                    'chr1\t24\t-\tchr2\t11\t+',
                    'chr1\t4\t+\tchr2\t32\t-',
                    'chr1\t44\t+\tchr1\t16\t-',
                    'chr3\t10\t-\tchr2\t12\t-',
                    'chr1\t20\t+\tchr3\t20\t-',
                    'chr1\t12\t-\tchr2\t26\t-',
                    'chr1\t14\t+\tchr1\t37\t+',
                    'chr1\t1\t+\tchr2\t10\t+'
                ))
   
    def tearDown(self):
        ''' Remove temporary files and directories '''
        for file in [self.chrFile, self.inBed1, self.inBed2, self.inFrag]:
            if os.path.isfile(file):
                os.remove(file)
        os.removedirs(self.dirName)

class TestBinFormation(InteractionTestCase):

    def test_bin_generation1(self):
        ''' Test unequal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,10,False))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([1,11,20,29,38], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([10,19,28,37,46], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1,2,3,4], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([1,10,18,26], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([9,17,25,33], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([5,6,7,8], dtype=np.uint32)))
       
   
    def test_bin_generation2(self):
        ''' Test unequal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,11,False))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([1,11,20,29,38], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([10,19,28,37,46], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1,2,3,4], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([1,12,23], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([11,22,33], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([5,6,7], dtype=np.uint32)))
   
    def test_bin_generation3(self):
        ''' Test unequal bin creation with large bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,40,False))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([1,24], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([23,46], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([1], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([33], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([2], dtype=np.uint32)))
   
    def test_bin_generation4(self):
        ''' Test equal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,10,True))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([4,14,24,34], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([13,23,33,43], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1,2,3], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([2,12,22], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([11,21,31], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([4,5,6], dtype=np.uint32)))
   
    def test_bin_generation5(self):
        ''' Test equal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,11,True))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([2,13,24,35], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([12,23,34,45], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1,2,3], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([1,12,23], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([11,22,33], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([4,5,6], dtype=np.uint32)))
   
    def test_bin_generation6(self):
        ''' Test unequal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,41,True))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([2,13,24,35], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([12,23,34,45], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1,2,3], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([1,12,23], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([11,22,33], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([4,5,6], dtype=np.uint32)))
     
    def test_bin_generation6(self):   
        ''' Test unequal bin creation with small bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,41,True))
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([3], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([43], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([], dtype=np.uint32)))
    
    def test_bin_generation7(self):   
        ''' Test bin creation from bed file '''
        genomeBin = interactionMatrix.genomeBin(self.inBed1)
        self.assertTrue(all(genomeBin.binDict['chr1']['start'] ==
            np.array([10,30], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['end'] ==
            np.array([20,40], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr1']['index'] ==
            np.array([0,1], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['start'] ==
            np.array([5], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['end'] ==
            np.array([15], dtype=np.uint32)))
        self.assertTrue(all(genomeBin.binDict['chr2']['index'] ==
            np.array([2], dtype=np.uint32)))
    
    def test_bin_generation8(self):
        ''' Test overlapping bins in bed file '''
        with self.assertRaises(IOError):
            interactionMatrix.genomeBin(self.inBed2)

class TestIndexFinder(InteractionTestCase):

    def test_index_finder1(self):
        ''' Test finding non contiguous bins '''
        genomeBin = interactionMatrix.genomeBin(self.inBed1)
        self.assertEqual(genomeBin.findBinIndex('chr1',10), 0)
        self.assertEqual(genomeBin.findBinIndex('chr1',9), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr1',30), 1)
        self.assertEqual(genomeBin.findBinIndex('chr1',29), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr1',20), 0)
        self.assertEqual(genomeBin.findBinIndex('chr1',21), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr1',40), 1)
        self.assertEqual(genomeBin.findBinIndex('chr1',41), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr3',10), 'nochr')
   
    def test_index_finder2(self):
        ''' Test finding contiguous unequal bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,10,False))
        self.assertEqual(genomeBin.findBinIndex('chr2',0), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr2',1), 5)
        self.assertEqual(genomeBin.findBinIndex('chr2',9), 5)
        self.assertEqual(genomeBin.findBinIndex('chr2',10), 6)
        self.assertEqual(genomeBin.findBinIndex('chr2',25), 7)
        self.assertEqual(genomeBin.findBinIndex('chr2',26), 8)
        self.assertEqual(genomeBin.findBinIndex('chr2',33), 8)
        self.assertEqual(genomeBin.findBinIndex('chr2',34), 'nobin')
   
    def test_index_finder3(self):
        ''' Test finding contiguous equal bins '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,10,True))
        self.assertEqual(genomeBin.findBinIndex('chr1',1), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr2',32), 'nobin')
        self.assertEqual(genomeBin.findBinIndex('chr2',31), 6)

class TestMatrixGeneration(InteractionTestCase):
   
    def test_matrix_generation(self):
        ''' Test creation of matrix '''
        genomeBin = interactionMatrix.genomeBin((self.chrFile,10,True))
        countMatrix, logArray = genomeBin.generateMatrix(
            self.inFrag, threads=4)
        self.assertTrue(np.array_equal(countMatrix,
            np.array([
                [2,0,0,0,0,0,1],
                [0,0,0,1,0,0,0],
                [0,0,0,0,2,0,0],
                [0,1,0,0,0,0,0],
                [0,0,2,0,0,0,0],
                [0,0,0,0,0,0,0],
                [1,0,0,0,0,0,0]
            ])))
        self.assertTrue(np.array_equal(logArray, np.array([10,2,3,5])))

suite = unittest.TestLoader().loadTestsFromTestCase(TestBinFormation)
unittest.TextTestRunner(verbosity=3).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestIndexFinder)
unittest.TextTestRunner(verbosity=3).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestMatrixGeneration)
unittest.TextTestRunner(verbosity=3).run(suite)
