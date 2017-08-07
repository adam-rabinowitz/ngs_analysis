from banyan import SortedSet, OverlappingIntervalsUpdator

class IntervalTree(object):
    
    def __init__(
            self
        ):
        # Create tree
        self.tree = SortedSet(updator=OverlappingIntervalsUpdator)
    
    def check_interval(
            self, start, end
        ):
        ''' Check intervals '''
        # Check arguments
        if not isinstance(start, int):
            raise TypeError('start must be an integer')
        if not isinstance(end, int):
            raise TypeError('end must be an integer')
        if end < start:
            raise ValueError('end is less than start')
    
    def add(
            self, start, end, data=None
        ):
        ''' Add intervals to tree '''
        # Check interval and add to tree
        self.check_interval(start, end)
        self.tree.add((start, end, data))
    
    def overlap(
            self, start, end
        ):
        ''' Find overlapping intervals '''
        # Check interval and find overlaps
        self.check_interval(start, end)
        width = end - start
        # Generate overlaps for intervals of 2 or more bases
        if width > 1:
            start += 1
            end -= 1
            overlaps = self.tree.overlap((start, end))
        # Generate and filter overlaps for shoter intervals
        else:
            overlaps = self.tree.overlap((start, end))
            overlaps = [x for x in overlaps if x[1] > start and x[0] < end]
        # Return data
        return(overlaps)
    
    def left(
            self, position, max_dist
        ):
        ''' Find leftermost intervals '''
        # Check arguments and find overlaps
        start = position - max_dist
        end = position
        self.check_interval(start, end)
        # Get overlaps and sort by end
        overlaps = self.tree.overlap((start, end))
        overlaps = sorted(overlaps, key=lambda x: (x[1], -x[0]))
        # Create output
        output = []
        while overlaps:
            iv = overlaps.pop()
            if iv[1] <= end:
                if output:
                    if iv[1] < output[0][1]:
                        break
                output.append(iv)
        # Return output
        return(output)
    
    def right(
            self, position, max_dist
        ):
        ''' Find rightermost intervals '''
        # Check arguments and find overlaps
        start = position
        end = position + max_dist
        self.check_interval(start, end)
        # Get overlaps and sort by end
        overlaps = self.tree.overlap((start, end))
        overlaps = sorted(overlaps, key=lambda x: (-x[0], -x[1]))
        # Create output
        output = []
        while overlaps:
            iv = overlaps.pop()
            if iv[0] >= start:
                if output:
                    if iv[0] > output[0][0]:
                        break
                output.append(iv)
        # Return output
        return(output)

# Perform unittesting
if __name__ == "__main__":
    
    import unittest
    
    class FindIntervalTest(unittest.TestCase):
        
        def setUp(self):
            # Create intervals for unittest
            self.int1  = (30, 50, {'name':'A'})
            self.int2  = (35, 55, {'name':'B'})
            self.int3  = (40, 60, {'name':'C'})
            self.int4  = (45, 65, {'name':'D'})
            self.int5  = (50, 70, {'name':'E'})
            self.int6  = (55, 75, {'name':'F'})
            self.int7  = (60, 80, {'name':'G'})
            self.int8  = (65, 85, {'name':'H'})
            self.int9  = (70, 90, {'name':'I'})
            self.int10 = (30, 30, {'name':'J'})
            self.int11 = (60, 60, {'name':'K'})
            self.int12 = (90, 90, {'name':'L'})
            # Build interval object for testing
            self.it = IntervalTree()
            for interval in (
                    self.int1, self.int2, self.int3, self.int4, self.int5,
                    self.int6, self.int7, self.int8, self.int9, self.int10,
                    self.int11, self.int12
                ):
                self.it.add(*interval)
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
                return(True)
            self.tup_equal = tup_equal
    
    class FindOverlapTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
            
        def test_wide_overlap_1(self):
            overlap = self.it.overlap(50, 70)
            expected = [self.int2, self.int3, self.int4, self.int5, self.int6,
                self.int11, self.int7, self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_2(self):
            overlap = self.it.overlap(49, 71)
            expected = [self.int1, self.int2, self.int3, self.int4, self.int5,
                self.int6, self.int11, self.int7, self.int8, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_3(self):
            overlap = self.it.overlap(30, 90)
            expected = [self.int1, self.int2, self.int3, self.int4, self.int5,
                self.int6, self.int11, self.int7, self.int8, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_wide_overlap_4(self):
            overlap = self.it.overlap(29, 91)
            expected = [self.int10, self.int1, self.int2, self.int3, self.int4,
                self.int5, self.int6, self.int11, self.int7, self.int8,
                self.int9, self.int12]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_1(self):
            overlap = self.it.overlap(60, 60)
            expected = [self.int4, self.int5, self.int6]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_2(self):
            overlap = self.it.overlap(60, 61)
            expected = [self.int4, self.int5, self.int6, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_3(self):
            overlap = self.it.overlap(59, 60)
            expected = [self.int3, self.int4, self.int5, self.int6]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_narrow_overlap_4(self):
            overlap = self.it.overlap(59, 61)
            expected = [self.int3, self.int4, self.int5, self.int6,
                self.int11, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
      
    class FindLeftTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
        
        def test_left_1(self):
            overlap = self.it.left(100, 10)
            expected = [self.int9, self.int12]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_2(self):
            overlap = self.it.left(100, 9)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_3(self):
            overlap = self.it.left(60, 10)
            expected = [self.int3, self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_4(self):
            overlap = self.it.left(61, 0)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_left_5(self):
            overlap = self.it.left(59, 10)
            expected = [self.int2]
            self.assertTrue(self.tup_equal(overlap, expected))
    
    class FindRightTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
        
        def test_right_1(self):
            overlap = self.it.right(20, 10)
            expected = [self.int10, self.int1]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_2(self):
            overlap = self.it.right(20, 9)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_3(self):
            overlap = self.it.right(60, 10)
            expected = [self.int11, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_4(self):
            overlap = self.it.right(59, 0)
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_right_5(self):
            overlap = self.it.right(61, 10)
            expected = [self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
    
    class CheckErrorTest(FindIntervalTest):
        
        def test_add_error_1(self):
            with self.assertRaises(TypeError):
                self.it.add(10, 'B')
            with self.assertRaises(TypeError):
                self.it.add('B', 10)
            with self.assertRaises(TypeError):
                self.it.add('B', 'C')
        
        def test_add_error_2(self):
            with self.assertRaises(ValueError):
                self.it.add(10, 9)
        
        def test_overlap_error_1(self):
            with self.assertRaises(TypeError):
                self.it.overlap(10, 'B')
            with self.assertRaises(TypeError):
                self.it.overlap('B', 10)
            with self.assertRaises(TypeError):
                self.it.overlap('B', 'C')
        
        def test_overlap_error_2(self):
            with self.assertRaises(ValueError):
                self.it.overlap(10, 9)
        
        def test_left_error_1(self):
            with self.assertRaises(TypeError):
                self.it.left(10, 'B')
            with self.assertRaises(TypeError):
                self.it.left('B', 10)
            with self.assertRaises(TypeError):
                self.it.left('B', 'C')
        
        def test_left_error_2(self):
            with self.assertRaises(ValueError):
                self.it.left(10, -1)
        
        def test_right_error_1(self):
            with self.assertRaises(TypeError):
                self.it.right(10, 'B')
            with self.assertRaises(TypeError):
                self.it.right('B', 10)
            with self.assertRaises(TypeError):
                self.it.right('B', 'C')
        
        def test_right_error_2(self):
            with self.assertRaises(ValueError):
                self.it.right(10, -1)
        
    suite = unittest.TestLoader().loadTestsFromTestCase(FindOverlapTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindLeftTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindRightTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(CheckErrorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
