# Import modules
import unittest
import os
import tempfile
import numpy
import collections
# Import personal modules
from chromosome_structure import fragend


class FragendTestCase(unittest.TestCase):
    
    def setUp(self):
        ''' Create required data for unittest '''
        # Create directories and file names
        self.dirName = tempfile.mkdtemp()
        self.fastaIn = self.dirName + '/test.fasta'
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
            ('chr4','+'):[],('chr4','-'):[]
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
        # Create files
        testIn = open(self.fastaIn, 'w')
        testIn.write('>chr1\n%s\n>chr2\n%s\n>chr3\n%s\n>chr4\n%s' %(
            self.chr1,
            self.chr2,
            self.chr3,
            self.chr4
        ))
        testIn.close()
    
    def tearDown(self):
        ''' Remove temporary files and directories '''
        if os.path.isfile(self.fastaIn):
            os.remove(self.fastaIn)
        os.removedirs(self.dirName)


class TestFragendAnalysis(FragendTestCase):
    ''' Class to test Fragend objects '''
    
    def test_creation(self):
        ''' Test creation of fragend object '''
        fe = fragend.Fragend(fasta = self.fastaIn, resite = self.resite)
        self.assertEqual(
            fe.frags,
            self.fragDict
        )
        self.assertEqual(
            fe.resite,
            'GATC'
        )
        self.assertEqual(
            fe.chrom,
            ['chr1','chr2','chr3','chr4']
        )
    
    def test_forward_strand_fragend(self):
        ''' Test identification of downstream fragends on + strand. '''
        fe = fragend.Fragend(fasta = self.fastaIn, resite = self.resite)
        self.assertEqual(
            fe.downstream([self.read1]),
            [('chr1',10,'+',8)]
        )
        self.assertEqual(
            fe.downstream([self.read2]),
            [('chr1',24,'+',21)]
        )
        self.assertEqual(
            fe.downstream([self.read3]),
            [('chr3',32,'+',8)]
        )
        self.assertEqual(
            fe.downstream([self.read4]),
            [('chr1',None,'+',None)]
        )
        self.assertEqual(
            fe.downstream([self.read5]),
            [('chr4',None,'+',None)]
        )
    
    def test_reverse_strand_fragend(self):
        ''' Test identification of downstream fragends on - strand. '''
        fe = fragend.Fragend(fasta = self.fastaIn, resite = self.resite)
        self.assertEqual(
            fe.downstream([self.read6]),
            [('chr2',20,'-',8)]
        )
        self.assertEqual(
            fe.downstream([self.read7]),
            [('chr2',16,'-',11)]
        )
        self.assertEqual(
            fe.downstream([self.read8]),
            [('chr3',3,'-',8)]
        )
        self.assertEqual(
            fe.downstream([self.read9]),
            [('chr2',None,'-',None)]
        )
        self.assertEqual(
            fe.downstream([self.read10]),
            [('chr4',None,'-',None)]
        )
    
    def test_errors_fragend(self):
        ''' Test error reporting. '''
        fe = fragend.Fragend(fasta = self.fastaIn, resite = self.resite)
        with self.assertRaises(ValueError):
            fe.downstream([('chr5','10','19','+')])
        with self.assertRaises(ValueError):
            fe.downstream([('chr1','19','10','-')])
    
    def test_fragend_pairs(self):
        ''' Test identification of acceptable fragend pairs. '''
        fe = fragend.Fragend(fasta = self.fastaIn, resite = self.resite)
        pairData = fe.fragend_pairs(
            pairIn = [
                self.read1 + self.read3,
                self.read1 + self.read2,
                self.read6 + self.read7,
                self.read1 + self.read4,
                self.read3 + self.read8,
                self.read9 + self.read6
            ],
            maxDistance = 20
        )
        self.assertEqual(
            pairData[0],
            [('chr1',10,'+','chr3',32,'+'),
             ('chr2',20,'-','chr2',16,'-'),
             ('chr3',32,'+','chr3',3,'-')]
        )
        self.assertEqual(
            pairData[1],
            collections.OrderedDict([
                ('total', 6),
                ('no fragend', 2),
                ('fragend too distant', 1),
                ('intrachromosomal', 2),
                ('interchromosomal', 1),
                ('fragend distance', [8,8,8,21,8,11,8,8]),
                ('ligation distance', [4,29])
            ])
        )


suite = unittest.TestLoader().loadTestsFromTestCase(TestFragendAnalysis)
unittest.TextTestRunner(verbosity=3).run(suite)
