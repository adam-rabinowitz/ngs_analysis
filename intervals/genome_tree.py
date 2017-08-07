import collections
from banyan_tree import IntervalTree

# Create new interval tree class
class CheckFreeTree(IntervalTree):
    
    def check_interval(
            self, start, end
        ):
        pass

# Create genome tree class
class GenomeTree(object):
    
    def __init__(self, ref_lengths):
        # Create variables to store data
        self.tree = collections.OrderedDict()
        self.length = collections.OrderedDict()
        self.strand = ('+', '-', '.')
        # Loop through chromosome and lengths
        for reference, length in ref_lengths:
            # Check input data
            if not isinstance(reference, str):
                raise TypeError('reference must be strings')
            if not isinstance(length, int):
                raise TypeError('length must be an integer')
            if length < 1:
                raise ValueError('length must be > 0')
            # Store data
            for s in self.strand:
                self.tree[(reference, s)] = CheckFreeTree()
            self.length[reference] = length
    
    def check_interval(
            self, reference, start, end, strand
        ):
        # Check reference
        if not isinstance(reference, str):
            raise TypeError('reference must be a string')
        try:
            length = self.length[reference]
        except KeyError:
            raise ValueError('unknown reference name')
        # Check start, end and reference
        if not isinstance(start, int):
            raise TypeError('start must be an integer')
        if not isinstance(end, int):
            raise TypeError('end must be an integer')
        if start < 0:
            raise ValueError('start cannot be negative')
        if end < start:
            raise ValueError('end is less than start')
        if end > length:
            raise ValueError('end is greater than reference length')
        # Check strand
        if len(strand) > len(set(strand)):
            raise ValueError('repeat strand arguments')
        if not isinstance(strand, str):
            raise TypeError('')
        for s in strand:
            if s not in self.strand:
                raise ValueError('unknown strand')
    
    def add(
            self, reference, start, end, strand='.', data=None
        ):
        # Check arguments
        self.check_interval(reference, start, end, strand)
        if len(strand) != 1:
            raise ValueError('add only accepts single strand argument')
        # Add data
        self.tree[(reference, strand)].add(start, end, data)
    
    def overlap(
            self, reference, start, end, strand='+-.', sort=True
        ):
        # Check arguments
        self.check_interval(reference, start, end, strand)
        # Extract overlaps
        output = []
        for s in strand:
            overlaps = self.tree[(reference, s)].overlap(start, end)
            overlaps = [(reference, x[0], x[1], s, x[2]) for x in overlaps]
            output.extend(overlaps)
        # Sort and return data
        if sort:
            output = sorted(output, key=lambda x: (x[1], x[2], x[3]))
        return(output)
    
    def left(
            self, reference, position, max_dist, strand='+-.', sort=True
        ):
        # Check arguments
        self.check_interval(reference, position, position, strand)
        if not isinstance(max_dist, int):
            raise TypeError('max_dist not integer')
        if max_dist < 0:
            raise ValueError('max_dist is negative')
        # Extract overlaps
        output = []
        for s in strand:
            # Find and process overlaps
            left = self.tree[(reference, s)].left(position, max_dist)
            left = [(reference, x[0], x[1], s, x[2]) for x in left]
            # Extract distance
            try:
                distance = position - left[0][2]
            except IndexError:
                continue
            # Process distance and store values
            if distance < max_dist:
                output = left
                max_dist = distance
            else:
                output.extend(left)
        # Sort and return data
        if sort:
            output = sorted(output, key=lambda x: x[1:4])
        return(output)
    
    def right(
            self, reference, position, max_dist, strand='+-.', sort=True
        ):
        # Check arguments
        self.check_interval(reference, position, position, strand)
        if not isinstance(max_dist, int):
            raise TypeError('max_dist not integer')
        if max_dist < 0:
            raise ValueError('max_dist is negative')
        # Extract overlaps
        output = []
        for s in strand:
            # Find and process overlaps
            right = self.tree[(reference, s)].right(position, max_dist)
            right = [(reference, x[0], x[1], s, x[2]) for x in right]
            # Extract distance
            try:
                distance = right[0][1] - position
            except IndexError:
                continue
            # Process distance and store values
            if distance < max_dist:
                output = right
                max_dist = distance
            else:
                output.extend(right)
        # Sort and return data
        if sort:
            output = sorted(output, key=lambda x: x[1:4])
        return(output)

# Perform unittesting
if __name__ == "__main__":
    
    import unittest
    
    class FindIntervalTest(unittest.TestCase):
        
        def setUp(self):
            # Create intervals for unittest
            self.int1  = ('chr1', 30, 50, '+', {'name':'A'})
            self.int2  = ('chr1', 35, 55, '-', {'name':'B'})
            self.int3  = ('chr1', 40, 60, '.', {'name':'C'})
            self.int4  = ('chr1', 45, 65, '+', {'name':'D'})
            self.int5  = ('chr1', 50, 70, '-', {'name':'E'})
            self.int6  = ('chr1', 55, 75, '.', {'name':'F'})
            self.int7  = ('chr1', 60, 80, '+', {'name':'G'})
            self.int8  = ('chr1', 65, 85, '-', {'name':'H'})
            self.int9  = ('chr1', 70, 90, '.', {'name':'I'})
            self.int10 = ('chr1', 30, 30, '.', {'name':'J'})
            self.int11 = ('chr1', 60, 60, '-', {'name':'K'})
            self.int12 = ('chr1', 90, 90, '+', {'name':'L'})
            # Build interval object for testing
            self.gt = GenomeTree([('chr1', 120), ('chr2', 100)])
            for interval in (
                    self.int1, self.int2, self.int3, self.int4, self.int5,
                    self.int6, self.int7, self.int8, self.int9, self.int10,
                    self.int11, self.int12
                ):
                self.gt.add(*interval)
            # Function to compare list of tuples
            def tup_equal(list1, list2):
                if len(list1) != len(list2):
                    return(False)
                for tup1, tup2 in zip(list1, list2):
                    if tup1[0] != tup2[0]:
                        return(False)
                    if tup1[1] != tup2[1]:
                        return(False)
                    if tup1[2] != tup2[2]:
                        return(False)
                    if tup1[3] != tup2[3]:
                        return(False)
                    if tup1[4] != tup2[4]:
                        return(False)
                return(True)
            self.tup_equal = tup_equal
    
    class FindOverlapTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
            
        def test_wide_overlap_1(self):
            overlap = self.gt.overlap('chr1', 50, 70)
            expected = [self.int2, self.int3, self.int4, self.int5, self.int6,
                self.int11, self.int7, self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_2(self):
            overlap = self.gt.overlap('chr1', 49, 71, '-.+')
            expected = [self.int1, self.int2, self.int3, self.int4, self.int5,
                self.int6, self.int11, self.int7, self.int8, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_3(self):
            overlap = self.gt.overlap('chr1', 49, 71, '-+')
            expected = [self.int1, self.int2, self.int4, self.int5, self.int11,
                self.int7, self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_4(self):
            overlap = self.gt.overlap('chr1', 49, 71, '.')
            expected = [self.int3, self.int6, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_1(self):
            overlap = self.gt.overlap('chr1', 60, 60, '+-.')
            expected = [self.int4, self.int5, self.int6]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_2(self):
            overlap = self.gt.overlap('chr1', 59, 61)
            expected = [self.int3, self.int4, self.int5, self.int6,
                self.int11, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_3(self):
            overlap = self.gt.overlap('chr1', 59, 61, '-.')
            expected = [self.int3, self.int5, self.int6, self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
    
        def test_narrow_overlap_4(self):
            overlap = self.gt.overlap('chr1', 59, 61, '+')
            expected = [self.int4, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
    
    class FindLeftTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
        
        def test_left_1(self):
            overlap = self.gt.left('chr1', 60, 10)
            expected = [self.int3, self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_2(self):
            overlap = self.gt.left('chr1', 61, 10, '-')
            expected = [self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_3(self):
            overlap = self.gt.left('chr1', 61, 10, '+')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_4(self):
            overlap = self.gt.left('chr1', 61, 0)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_5(self):
            overlap = self.gt.left('chr1', 59, 10)
            expected = [self.int2]
            self.assertTrue(self.tup_equal(overlap, expected))
    
    class FindRightTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
        
        def test_right_1(self):
            overlap = self.gt.right('chr1', 60, 10)
            expected = [self.int11, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_2(self):
            overlap = self.gt.right('chr1', 59, 10, '+')
            expected = [self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_3(self):
            overlap = self.gt.right('chr1', 59, 0, '.')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_4(self):
            overlap = self.gt.right('chr1', 59, 0)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_5(self):
            overlap = self.gt.right('chr1', 61, 10)
            expected = [self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
    
    class CheckErrorTest(FindIntervalTest):
        
        def test_add_type_error(self):
            with self.assertRaises(TypeError):
                self.gt.add(0, 10, 20, '+-.')
            with self.assertRaises(TypeError):
                self.gt.add('chr1', '10', 20, '+-.')
            with self.assertRaises(TypeError):
                self.gt.add('chr1', 10, '20', '+-.')
            with self.assertRaises(TypeError):
                self.gt.add('chr1', 10, 20, 0)
        
        def test_add_value_error(self):
            with self.assertRaises(ValueError):
                self.gt.add('chr3', 10, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.add('chr1', -1, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.add('chr1', 10, 121, '+-.')
            with self.assertRaises(ValueError):
                self.gt.add('chr1', 11, 10, '+-.')
            with self.assertRaises(ValueError):
                self.gt.add('chr1', 10, 20, '*')
            with self.assertRaises(ValueError):
                self.gt.add('chr1', 10, 20, '++')
        
        def test_overlap_error(self):
            with self.assertRaises(ValueError):
                self.gt.overlap('chr3', 10, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.overlap('chr1', -1, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.overlap('chr1', 10, 121, '+-.')
            with self.assertRaises(ValueError):
                self.gt.overlap('chr1', 11, 10, '+-.')
            with self.assertRaises(ValueError):
                self.gt.overlap('chr1', 10, 20, '*')
            with self.assertRaises(ValueError):
                self.gt.overlap('chr1', 10, 20, '++')
        
        def test_left_error(self):
            with self.assertRaises(ValueError):
                self.gt.left('chr3', 10, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.left('chr1', -1, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.left('chr1', 121, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.left('chr1', 10, -1, '+')
            with self.assertRaises(ValueError):
                self.gt.left('chr1', 10, 20, '*')
            with self.assertRaises(ValueError):
                self.gt.left('chr1', 10, 20, '++')
        
        def test_right_error(self):
            with self.assertRaises(ValueError):
                self.gt.right('chr3', 10, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.right('chr1', -1, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.right('chr1', 121, 20, '+-.')
            with self.assertRaises(ValueError):
                self.gt.right('chr1', 10, -1, '+-.')
            with self.assertRaises(ValueError):
                self.gt.right('chr1', 10, 20, '*')
            with self.assertRaises(ValueError):
                self.gt.right('chr1', 10, 20, '++')
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindOverlapTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindLeftTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindRightTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(CheckErrorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
