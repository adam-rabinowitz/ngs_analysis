import unittest
import tempfile
import collections
import os
import multiprocessing
from ngs_analysis.structure import alignedPair

class PairTestCase(unittest.TestCase):
    
    def setUp(self):
        ''' Create temporary directory and example read pairs '''
        # Make temporary file
        self.dirName = tempfile.mkdtemp()
        self.testPair = self.dirName + 'test.pair'
        # Create read pairs
        self.pair1 = ('chr1',1,40,'+','chr1',1960,2000,'-')
        self.pair2 = ('chr1',1,40,'+','chr1',1959,2001,'-')
        self.pair3 = ('chr1',1,40,'+','chr2',1959,1999,'-')
        self.pair4 = ('chr1',1,40,'+','chr1',1959,1999,'+')
        self.pair5 = ('chr1',100,140,'-','chr1',100,140,'+')
        self.pair6 = ('chr1',100,140,'-','chr1',90,130,'+')
        self.pair7 = ('chr1',100,140,'-','chr1',90,141,'+')
        self.pair8 = ('chr1',99,140,'-','chr1',100,130,'+')
        # Create pair list
        self.pairList = ([self.pair2] + [self.pair3] * 2 + [self.pair4] * 3 +
            [self.pair5] + [self.pair6] * 2)
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        if os.path.isfile(self.testPair):
            os.remove(self.testPair)
        os.removedirs(self.dirName)
    
    def readFile(self):
        with open(self.testPair) as f:
            data = f.readlines()
        output = [d.strip().split('\t') for d in data]
        return(output)
    
    def processPair(self, rmDup, rmConcord, maxSize):
        pipes = multiprocessing.Pipe(True)
        process = multiprocessing.Process(
            target = alignedPair.processPairs,
            args = (pipes[0], self.testPair, rmDup, rmConcord, maxSize)
        )
        process.start()
        pipes[0].close()
        for pair in self.pairList:
            pipes[1].send(pair)
        pipes[1].send(None)
        metrics = pipes[1].recv()
        pipes[1].close()
        process.join()
        return(metrics)


class TestPairProcessing(PairTestCase):
    
    def test_find_concordant(self):
        ''' Testing identification of concordant read pairs '''
        # Check proper pair
        self.assertTrue(alignedPair.concordant(self.pair1,2000))
        # Check pair that is too big
        self.assertFalse(alignedPair.concordant(self.pair2,2000))
        # Check pair on different chromosome
        self.assertFalse(alignedPair.concordant(self.pair3,2000))
        # Check pair on same strand
        self.assertFalse(alignedPair.concordant(self.pair4,2000))
        # Check overlapping proper pairs
        self.assertTrue(alignedPair.concordant(self.pair5,2000))
        self.assertTrue(alignedPair.concordant(self.pair6,2000))
        # Check when read pairs extend beyond each other
        self.assertFalse(alignedPair.concordant(self.pair7,2000))
        self.assertFalse(alignedPair.concordant(self.pair8,2000))
    
    def test_process_concord_duplication(self):
        ''' Test correct processing of concordant and duplicated reads '''
        # Check processing with concordant and duplicates removed
        pairMetrics = self.processPair(rmDup = True, rmConcord = True,
            maxSize = 2000)
        self.assertEqual(self.readFile(), [map(str,self.pair2),
            map(str,self.pair3), map(str,self.pair4)])
        self.assertEqual(pairMetrics, collections.defaultdict(int, {
            'total':9, 'unique':5, 'duplicate':4, 'concord':3, 'concorduni':2,
            'discord':6, 'discorduni':3}))
        # Check processing with duplicates removed
        pairMetrics = self.processPair(rmDup = True, rmConcord = False,
            maxSize = 2000)
        self.assertEqual(self.readFile(), [map(str,self.pair2),
            map(str,self.pair3), map(str,self.pair4), map(str,self.pair5),
            map(str,self.pair6)])
        self.assertEqual(pairMetrics, collections.defaultdict(int, {
            'total':9, 'unique':5, 'duplicate':4, 'concord':3, 'concorduni':2,
            'discord':6, 'discorduni':3}))
        # Check processing with concordant removed
        pairMetrics = self.processPair(rmDup = False, rmConcord = True,
            maxSize = 2000)
        self.assertEqual(self.readFile(), [map(str,self.pair2)] +
            [map(str,self.pair3)] * 2 + [map(str,self.pair4)] * 3)
        self.assertEqual(pairMetrics, collections.defaultdict(int, {
            'total':9, 'unique':5, 'duplicate':4, 'concord':3, 'concorduni':2,
            'discord':6, 'discorduni':3}))
        # Check processing with nothing removed
        pairMetrics = self.processPair(rmDup = False, rmConcord = False,
            maxSize = 2000)
        self.assertEqual(self.readFile(), [map(str,self.pair2)] +
            [map(str,self.pair3)] * 2 + [map(str,self.pair4)] * 3 +
            [map(str,self.pair5)] + [map(str,self.pair6)] * 2)
        self.assertEqual(pairMetrics, collections.defaultdict(int, {
            'total':9, 'unique':5, 'duplicate':4, 'concord':3, 'concorduni':2,
            'discord':6, 'discorduni':3}))


suite = unittest.TestLoader().loadTestsFromTestCase(TestPairProcessing)
unittest.TextTestRunner(verbosity=3).run(suite)

