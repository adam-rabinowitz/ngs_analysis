import unittest
import pysam
from ngs_python.bam import pysamfunc

def compareDict(d1,d2):
    # Find combined keys
    keySet = set()
    keySet.update(d1.keys())
    keySet.update(d2.keys())
    # Loop through keys and find differing elements
    for key in keySet:
        try:
            if d1[key] != d2[key]:
                print 'Values different for %s' %(key)
                return(False)
        except KeyError:
            print 'Absent values for %s' %(key)
            return(False)
    return(True)

class baseCallTest(unittest.TestCase):
    
    def test_match(self):
        # Create read
        read = pysam.AlignedSegment()
        read.query_sequence = 'AGCTAGTATGTA'
        read.query_qualities = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
        read.reference_start = 10
        read.cigartuples = [(0,12)]
        # Create expected result
        expDict = {10:('A',20), 11:('G',21), 12:('C',22), 13:('T',23),
            14:('A',24), 15:('G',25), 16:('T',26), 17:('A',27), 18:('T',28),
            19:('G',29), 20:('T',30), 21:('A',31)}
        # Create actual result and compare
        baseDict = pysamfunc.baseCalls(read)
        self.assertTrue(baseDict == expDict)
    
    def test_insertion(self):
        # Create read
        read = pysam.AlignedSegment()
        read.query_sequence = 'AGCTAGTATGTA'
        read.query_qualities = [20, 21, 22, 23, 24, 27, 26, 27, 28, 29, 30, 31]
        read.reference_start = 10
        read.cigartuples = [(0,3),(1,2),(0,2),(1,2),(0,3)]
        # Create expected result
        expDict = {10:('A',20), 11:('G',21), 12:('CTA',23), 13:('G',27),
            14:('TAT',27), 15:('G',29), 16:('T',30), 17:('A',31)}
        # Create actual result and compare
        baseDict = pysamfunc.baseCalls(read)
        self.assertTrue(baseDict == expDict)
    
    def test_deletion(self):
        # Create read
        read = pysam.AlignedSegment()
        read.query_sequence = 'AGCTAGTATGTA'
        read.query_qualities = [20, 21, 22, 28, 24, 25, 26, 27, 32, 29, 30, 31]
        read.reference_start = 10
        read.cigartuples = [(0,3),(2,2),(0,5),(2,2),(0,4)]
        # Create expected result
        expDict = {10:('A',20), 11:('G',21), 12:('C',22), 13:('-',25),
            14:('-',25),15:('T',28), 16:('A',24), 17:('G',25), 18:('T',26),
            19:('A',27), 20:('-',29), 21:('-',29), 22:('T',32), 23:('G',29),
            24:('T',30), 25:('A',31)}
        # Create actual result and compare
        baseDict = pysamfunc.baseCalls(read)
        self.assertTrue(baseDict == expDict)
    
    def test_deletion_group(self):
        # Create read
        read = pysam.AlignedSegment()
        read.query_sequence = 'AGCTAGTATGTA'
        read.query_qualities = [20, 21, 22, 28, 24, 25, 26, 27, 32, 29, 30, 31]
        read.reference_start = 10
        read.cigartuples = [(0,3),(2,2),(0,5),(2,2),(0,4)]
        # Create expected result
        expDict = {10:('A',20), 11:('G',21), 12:('C',22), 13:('--',25),
            15:('T',28), 16:('A',24), 17:('G',25), 18:('T',26), 19:('A',27),
            20:('--',29), 22:('T',32), 23:('G',29), 24:('T',30), 25:('A',31)}
        # Create actual result and compare
        baseDict = pysamfunc.baseCalls(read, groupdel = True)
        self.assertTrue(baseDict == expDict)

suite = unittest.TestLoader().loadTestsFromTestCase(baseCallTest)
unittest.TextTestRunner(verbosity=2).run(suite)

class startIntervalTest(unittest.TestCase):
    
    def test_forward_interval(self):
        # Create read
        read = pysam.AlignedSegment()
        read.reference_start = 100
        read.query_sequence = 'AGCTAGCTAGCT'
        read.cigartuples = [(0,12)]
        read.is_reverse = False
        read.reference_id = 0
        chrDict = {0:('chr1',200)}
        # Create expected result and compare
        intervalSize = 10
        expTuple = ('chr1', 95, 105)
        actTuple = pysamfunc.startInterval(read, intervalSize, chrDict)
        self.assertTrue(actTuple == expTuple)
        # Create expected result and compare
        intervalSize = 50
        expTuple = ('chr1', 75, 125)
        actTuple = pysamfunc.startInterval(read, intervalSize, chrDict)
        self.assertTrue(actTuple == expTuple)
    
    def test_reverse_interval(self):
        # Create read
        read = pysam.AlignedSegment()
        read.reference_start = 100
        read.query_sequence = 'AGCTAGCTAGCT'
        read.cigartuples = [(0,12)]
        read.is_reverse = True
        read.reference_id = 0
        chrDict = {0:('chr1',200)}
        # Create expected result and compare
        intervalSize = 10
        expTuple = ('chr1', 107, 117)
        actTuple = pysamfunc.startInterval(read, intervalSize, chrDict)
        self.assertTrue(actTuple == expTuple)
        # Create expected result and compare
        intervalSize = 50
        expTuple = ('chr1', 87, 137)
        actTuple = pysamfunc.startInterval(read, intervalSize, chrDict)
        self.assertTrue(actTuple == expTuple)

    def test_interval_trim(self):
        # Create read
        read = pysam.AlignedSegment()
        read.reference_start = 10
        read.query_sequence = 'AGCTAGCTAG'
        read.cigartuples = [(0,10)]
        read.is_reverse = False
        read.reference_id = 0
        chrDict = {0:('chr1',20)}
        intervalSize = 30
        # Create expected result and compare
        expTuple = ('chr1', 0, 20)
        actTuple = pysamfunc.startInterval(read, intervalSize, chrDict)
        self.assertTrue(actTuple == expTuple)

suite = unittest.TestLoader().loadTestsFromTestCase(startIntervalTest)
unittest.TextTestRunner(verbosity=2).run(suite)
