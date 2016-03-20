import unittest
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
        baseDict = pysamfunc.baseCalls(
            ['A','G','C','T','A','G','T','A','T','G','T','A'],
            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
            [(0,12)]
        )
        expDict = {10:('A',20), 11:('G',21), 12:('C',22), 13:('T',23),
            14:('A',24), 15:('G',25), 16:('T',26), 17:('A',27), 18:('T',28),
            19:('G',29), 20:('T',30), 21:('A',31)}
        self.assertTrue(baseDict == expDict)
    
    def test_insertion(self):
        baseDict = pysamfunc.baseCalls(
            ['A','G','C','T','A','G','T','A','T','G','T','A'],
            [20, 21, 22, 23, 24, 27, 26, 27, 28, 29, 30, 31],
            [10, 11, 12, 13, 14, 15, 16, 17],
            [(0,3),(1,2),(0,2),(1,2),(0,3)]
        )
        expDict = {10:('A',20), 11:('G',21), 12:('CTA',23), 13:('G',27),
            14:('TAT',27), 15:('G',29), 16:('T',30), 17:('A',31)}
        self.assertTrue(baseDict == expDict)
    
    def test_deletion(self):
        baseDict = pysamfunc.baseCalls(
            ['A','G','C','T','A','G','T','A','T','G','T','A'],
            [20, 21, 22, 28, 24, 25, 26, 27, 33, 29, 30, 31],
            [10, 11, 12, 15, 16, 17, 18, 19, 22, 23, 24, 25],
            [(0,3),(2,2),(0,5),(2,2),(0,4)],
        )
        expDict = {10:('A',20), 11:('G',21), 12:('C',22), 13:('-',25),
            14:('-',25),15:('T',28), 16:('A',24), 17:('G',25), 18:('T',26),
            19:('A',27), 20:('-',30), 21:('-',30), 22:('T',33), 23:('G',29),
            24:('T',30), 25:('A',31)}
        self.assertTrue(baseDict == expDict)

suite = unittest.TestLoader().loadTestsFromTestCase(baseCallTest)
unittest.TextTestRunner(verbosity=2).run(suite)
