# Import modules
import unittest
import os
import tempfile
import numpy
import collections
# Import personal modules
from ngs_analysis.structure import fragendPair

class FragendTestCase(unittest.TestCase):
    
    def setUp(self):
        ''' Create required data for unittest '''
        # Create directories and file names
        self.dirName = tempfile.mkdtemp()
        self.fastaIn = self.dirName + '/test.fasta'
        self.pairIn = self.dirName + '/test.pair'
        self.fragendOut = self.dirName + '/test.fragend'
        # Create data
        self.chr1 = 'AAAAAAAAGATCAAAAAAAAAAGATCAAAAAAAA'
        self.chr2 = 'AAAAAAAAAAAAAGATCgatcAAAAAAAAAAAAA'
        self.chr3 = 'GatcAAAAAAAAAAAAAAAAAAAAAAAAAAGATC'
        self.chr4 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        self.resite = 'GatC'
        self.fragDict = {
            ('chr1','+'):[12,26],('chr1','-'):[9,23],
            ('chr2','+'):[17,21],('chr2','-'):[14,18],
            ('chr3','+'):[4,34],('chr3','-'):[1,31],
            ('chr4','+'):[],('chr4','-'):[],
            'resite':'GATC'
        }
        self.read1 = ['chr1','3','12','+']
        self.read2 = ['chr1','4','13','+']
        self.read3 = ['chr3','25','34','+']
        self.read4 = ['chr1','18','27','+']
        self.read5 = ['chr4','10','19','+']
        self.read6 = ['chr2','18','27','-']
        self.read7 = ['chr2','17','26','-']
        self.read8 = ['chr3','1','10','-']
        self.read9 = ['chr2','13','22','-']
        self.read10 = ['chr4','10','19','-']
        # Create fasta file
        fastaIn = open(self.fastaIn, 'w')
        fastaIn.write('>chr1\n%s\n>chr2\n%s\n>chr3\n%s\n>chr4\n%s' %(
            self.chr1,
            self.chr2,
            self.chr3,
            self.chr4
        ))
        fastaIn.close()
        # Create pair file
        pairIn = open(self.pairIn, 'w')
        pairIn.write('%s\n%s\n%s\n%s\n%s\n%s\n' %(
            '\t'.join(self.read1 + self.read3),
            '\t'.join(self.read1 + self.read2),
            '\t'.join(self.read6 + self.read7),
            '\t'.join(self.read1 + self.read4),
            '\t'.join(self.read3 + self.read8),
            '\t'.join(self.read9 + self.read6)
        ))
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        for file in [self.fastaIn, self.pairIn, self.fragendOut]:
            if os.path.isfile(file):
                os.remove(file)
        os.removedirs(self.dirName)


class TestFragendAnalysis(FragendTestCase):
    ''' Class to test Fragend objects '''
    
    def test_creation(self):
        ''' Test creation of fragend object '''
        fd = fragendPair.findFragendSites(fasta = self.fastaIn,
            resite = self.resite)
        self.assertEqual(fd, self.fragDict)
    
    def test_forward_strand_fragend(self):
        ''' Test identification of downstream fragends on + strand. '''
        self.assertEqual(
            fragendPair.downstream([self.read1], self.fragDict),
            [('chr1',10,'+',8)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read2], self.fragDict),
            [('chr1',24,'+',21)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read3], self.fragDict),
            [('chr3',32,'+',8)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read4], self.fragDict),
            [None]
        )
        self.assertEqual(
            fragendPair.downstream([self.read5], self.fragDict),
            [None]
        )
    
    def test_reverse_strand_fragend(self):
        ''' Test identification of downstream fragends on - strand. '''
        self.assertEqual(
            fragendPair.downstream([self.read6], self.fragDict),
            [('chr2',20,'-',8)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read7], self.fragDict),
            [('chr2',16,'-',11)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read8], self.fragDict),
            [('chr3',3,'-',8)]
        )
        self.assertEqual(
            fragendPair.downstream([self.read9], self.fragDict),
            [None]
        )
        self.assertEqual(
            fragendPair.downstream([self.read10], self.fragDict),
            [None]
        )
    
    def test_errors_fragend(self):
        ''' Test error reporting. '''
        with self.assertRaises(IOError):
            fragendPair.downstream([('chr5','10','19','+')], self.fragDict)
        with self.assertRaises(IOError):
            fragendPair.downstream([('chr1','19','10','-')], self.fragDict)
    
    def test_fragend_pairs(self):
        ''' Test identification of acceptable fragend pairs. '''
        pairData = fragendPair.fragendPairs(pairIn = self.pairIn,
            fasta = self.fastaIn, resite = 'gATc', maxDistance = 20,
            fragendOut = self.fragendOut)
        output = []
        with open(self.fragendOut) as fileIn:
            for line in fileIn:
                output.append(line.strip().split('\t'))
        self.assertEqual(
            output,
            [['chr1','10','+','chr3','32','+'],
             ['chr2','20','-','chr2','16','-'],
             ['chr3','32','+','chr3','3','-']]
        )
        self.assertEqual(
            pairData,
            collections.defaultdict(int, {'total':6, 'none':2, 'distant':1,
                'intrachromosomal':2, 'interchromosomal':1,
                'fragDist':[8,8,8,21,8,11,8,8], 'ligDist':[4,29]})
        )


suite = unittest.TestLoader().loadTestsFromTestCase(TestFragendAnalysis)
unittest.TextTestRunner(verbosity=3).run(suite)
