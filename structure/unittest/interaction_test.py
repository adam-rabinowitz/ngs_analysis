import unittest
import collections
import numpy as np
from chromosome_structure import interaction

class InteractionTestCase(unittest.TestCase):
    
    def setUp(self):
        ''' Create data for unittest.'''
        self.chrData = [['chr1',53],['chr2',47]]


class TestInteractionAnalysis(InteractionTestCase):

    def test_small_unqual_bins(self):
        ''' Test unequal bin creation with small bins '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 10
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['start'],
                np.array([1,10,19,28,37,46])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['end'],
                np.array([9,18,27,36,45,53])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['index'],
                np.array([0,1,2,3,4,5])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['start'],
                [1,11,21,30,39],
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['end'],
                np.array([10,20,29,38,47])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['index'],
                np.array([6,7,8,9,10])
            )
        )
        self.assertEqual(
            inter.binList,
            [('chr1',1,9),('chr1',10,18),('chr1',19,27),('chr1',28,36),
             ('chr1',37,45),('chr1',46,53),('chr2',1,10),('chr2',11,20),
             ('chr2',21,29),('chr2',30,38),('chr2',39,47)]
        )
        self.assertEqual(
            inter.binNames,
            ['chr1:1-9','chr1:10-18','chr1:19-27','chr1:28-36','chr1:37-45',
             'chr1:46-53','chr2:1-10','chr2:11-20','chr2:21-29','chr2:30-38',
             'chr2:39-47']
        )
    
    def test_large_unequal_bins(self):
        ''' Test unequal bin creation with large bins '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 100
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['start'],
                np.array([1])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['end'],
                np.array([53])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['index'],
                np.array([0])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['start'],
                np.array([1])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['end'],
                np.array([47])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['index'],
                np.array([1])
            )
        )
        self.assertEqual(
            inter.binList,
            [('chr1',1,53),('chr2',1,47)]
        )
        self.assertEqual(
            inter.binNames,
            ['chr1:1-53','chr2:1-47']
        )
    
    def test_small_equal_bins(self):
        ''' Test small equal bin formation '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 10,
            binEqual = True
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['start'],
                np.array([2,12,22,32,42]),
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['end'],
                np.array([11,21,31,41,51])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['index'],
                np.array([0,1,2,3,4])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['start'],
                np.array([4,14,24,34])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['end'],
                np.array([13,23,33,43])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['index'],
                np.array([5,6,7,8])
            )
        )
        self.assertEqual(
            inter.binList,
            [('chr1',2,11),('chr1',12,21),('chr1',22,31),('chr1',32,41),
             ('chr1',42,51),('chr2',4,13),('chr2',14,23),('chr2',24,33),
             ('chr2',34,43)]
        )
        self.assertEqual(
            inter.binNames,
            ['chr1:2-11','chr1:12-21','chr1:22-31','chr1:32-41','chr1:42-51',
             'chr2:4-13','chr2:14-23','chr2:24-33','chr2:34-43']
        )
    
    def test_large_equal_bins(self):
        ''' Test large equal bin formation. '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 100,
            binEqual = True
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['start'],
                np.array([])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['end'],
                np.array([])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr1']['index'],
                np.array([])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['start'],
                np.array([])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['end'],
                np.array([])
            )
        )
        self.assertTrue(
            np.array_equal(
                inter.chrBin['chr2']['index'],
                np.array([])
            )
        )
        self.assertEqual(
            inter.binList,
            []
        )
        self.assertEqual(
            inter.binNames,
            []
        )

    def test_find_bin_index(self):
        ''' Test small equal bin formation '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 10,
            binEqual = True
        )
        self.assertEqual(
            inter.index(
                [('chr1',2),('chr2',43),('chr1',1),('chr2',44),('chr1',32),
                 ('chr2',33),('chr1',31),('chr2',14),('chr2',4),('chr1',51)
                 ,('chr2',3),('chr1',52),('chr1',17)]
            ),
            [0,8,'','',3,7,2,6,5,4,'','',1]
        )

    def test_find_nobin_index(self):
        ''' Test finding of bins when chromosome has no bins '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 100,
            binEqual = True
        )
        self.assertEqual(
            inter.index([('chr1',25),('chr2',25)]),
            ['','']
        )


    def test_generate_count_array(self):
        ''' Test creation of a count array '''
        inter = interaction.Interaction(
            chrData = self.chrData,
            maxBinSize = 10,
            binEqual = True
        )
        array, log = inter.interactions(
            fragIn = [
                ('chr1',11,'+','chr2',43,'-'),
                ('chr1',2,'-','chr2',34,'+'),
                ('chr2',14,'+','chr2',23,'-'),
                ('chr1',26,'-','chr2',14,'+'),
                ('chr3',12,'+','chr1',15,'+'),
                ('chr1',34,'+','chr3',12,'+'),
                ('chr1',1,'-','chr2',23,'-'),
                ('chr1',11,'-','chr2',2,'+')
            ]
        )
        self.assertTrue(
            np.array_equal(
                array,
                np.array([
                    [0,0,0,0,0,0,0,0,2],
                    [0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0],
                    [0,0,1,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,0],
                    [2,0,0,0,0,0,0,0,0]
                ])
            )
        )
        self.assertEqual(
            log,
            collections.OrderedDict({
                'total':8,
                'no chromosome':2,
                'no bin':2,
                'accepted':4
            })
        )


suite = unittest.TestLoader().loadTestsFromTestCase(TestInteractionAnalysis)
unittest.TextTestRunner(verbosity=3).run(suite)
