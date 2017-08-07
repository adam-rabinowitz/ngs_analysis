from bx.intervals import Interval, Intersecter
import collections

class FindInterval:
    
    inttup = collections.namedtuple(
        'inttup', ['chrom', 'start', 'end', 'strand', 'data'])
    
    def __init__(
            self, chroms = None
        ):
        ''' Function to initialise object.
        
        Args:
        chroms (dict) - Optional dictionary of chromosome lengths.
        
        '''
        # Create and store key variables
        self.intervals = collections.defaultdict(Intersecter)
        if not chroms is None:
            if not isinstance(chroms, dict):
                raise TypeError('supplied chroms must be a dictionary')
        self.chroms = chroms
    
    def check_interval(
            self, chrom, start, end, strand
        ):
        ''' Function to check interval arguments.
        
        Args:
            chrom (str)- Name of chromosome.
            start (int)- Start of interval.
            end (int)- End of interval.
            strand (str)- Must be one of +/-/./*
            data (dict)- Dictionary of interval decriptors.
        
        Raises:
            TypeError - If argument is of wrong type.
            ValueError - If argument of wrong value.
        
        '''
        # Check all arguments
        if not isinstance(chrom, str):
            raise TypeError('chrom not a string')
        if not isinstance(start, int):
            raise TypeError('start not an integer')
        if not isinstance(end, int):
            raise TypeError('end not an integer')
        if start < 0:
            raise ValueError('start is negative')
        if end < start:
            raise ValueError('end is less than start')
        if not self.chroms is None:
            if end > self.chroms[chrom]:
                raise ValueError('end is larger than chromosome length')
        if not isinstance(strand, str):
            raise ValueError('strand is not a string')
        for s in strand:
            if strand.count(s) > 1:
                raise ValueError('strand value repeated')
            if s not in ('+-.'):
                print(s)
                raise ValueError('strand value not recognised')
    
    def check_position(
            self, chrom, position, strand
        ):
        ''' Function to check interval arguments.
        
        Args:
            chrom (str)- Name of chromosome.
            position (int)- Position on chromosome.
            strand (str)- Must be one of +/-/.
        
        Raises:
            TypeError - If argument is of wrong type.
            ValueError - If argument of wrong value.
        
        '''
        # Check all arguments
        if not isinstance(chrom, str):
            raise TypeError('chrom not a string')
        if not isinstance(position, int):
            raise TypeError('start not an integer')
        if position < 0:
            raise ValueError('pos is negative')
        if not self.chroms is None:
            if position > self.chroms[chrom]:
                raise ValueError('position is larger than chromosome length')
        if not isinstance(strand, str):
            raise TypeError('strand is not a string')
        for s in strand:
            if strand.count(s) > 1:
                raise ValueError('strand value repeated')
            if s not in ('+-.'):
                raise ValueError('strand value not recognised')
    
    def add_interval(
            self, chrom, start, end, strand, data={}
        ):
        ''' A function to add intervals to interval tree.
        
        Args:
            chrom (str)- Name of chromosome.
            start (int)- Start of interval.
            end (int)- End of interval.
            strand (str)- Must be one of +/-/.
            data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        self.check_interval(chrom, start, end, strand)
        if len(strand) > 1:
            raise ValueError('Ambiguous strand value supplied')
        if not isinstance(data, dict):
            raise TypeError('supplied data variable must be a dictionary')
        # Create interval and add to tree
        interval = Interval(start=start, end=end, value=data)
        self.intervals[(chrom, strand)].add_interval(interval)
    
    def find_overlap(
            self, chrom, start, end, strand='*'
        ):
        ''' A function to return data associated with overlapping intervals
        
        Args:
            chrom (str)- Name of chromosome.
            start (int)- Start of interval.
            end (int)- End of interval.
            strand (str)- Must be one of +/-/./*
        
        Returns:
            output - A list of named tuples containing five elements:
                chrom (str)- Name of chromosome.
                start (int)- Start of interval.
                end (int)- End of interval.
                strand (str)- Interval strand.
                data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        if strand == '*':
            strand = '+-.'
        self.check_interval(chrom, start, end, strand)
        # Find overlaps
        output = []
        for s in strand:
            for overlap in self.intervals[(chrom, s)].find(start, end):
                overtup = self.inttup(
                    chrom, overlap.start, overlap.end, s, overlap.value)
                output.append(overtup)
        # Sort overlaps and return
        output = sorted(output, key=lambda x: x[1:4])
        return(output)
         
    def find_start(
            self, chrom, position, strand = '*'
        ):
        ''' A function to return intervals with specified start
        
        Args:
            chrom (str)- Name of chromosome.
            position (int)- Position on chromosome.
            strand (str)- Must be a combination of '+/-/.' or '*'.
        
        Returns:
            output - A list of named tuples containing five elements:
                chrom (str)- Name of chromosome.
                start (int)- Start of interval.
                end (int)- End of interval.
                strand (str)- Interval strand.
                data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        if strand == '*':
            strand = '+-.'
        self.check_position(chrom, position, strand)
        # Find all intervals with end
        output = []
        for s in strand:
            for overlap in self.intervals[(chrom, s)].find(
                position - 1, position + 1):
                if overlap.start != position:
                    continue
                overtup = self.inttup(
                    chrom, overlap.start, overlap.end, s, overlap.value)
                output.append(overtup)
        # Sort overlaps and return
        output = sorted(output, key=lambda x: x[1:4])
        return(output)
    
    def find_end(
            self, chrom, position, strand = '*'
        ):
        ''' A function to return intervals with specified end
        
        Args:
            chrom (str)- Name of chromosome.
            position (int)- Position on chromosome.
            strand (str)- Must be a combination of '+/-/.' or '*'.
        
        Returns:
            output - A list of named tuples containing five elements:
                chrom (str)- Name of chromosome.
                start (int)- Start of interval.
                end (int)- End of interval.
                strand (str)- Interval strand.
                data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        if strand == '*':
            strand = '+-.'
        self.check_position(chrom, position, strand)
        # Find all intervals with end
        output = []
        for s in strand:
            for overlap in self.intervals[(chrom, s)].find(
                position - 1, position + 1):
                if overlap.end != position:
                    continue
                overtup = self.inttup(
                    chrom, overlap.start, overlap.end, s, overlap.value)
                output.append(overtup)
        # Sort overlaps and return
        output = sorted(output, key=lambda x: x[1:4])
        return(output)
    
    def find_left(
            self, chrom, position, strand = '*', dist=1000000000
        ):
        ''' A function to return data with leftmost intervals
        
        Args:
            chrom (str)- Name of chromosome.
            position (int)- Position on chromosome.
            strand (str)- Must be a combination of '+/-/.' or '*'.
            dist (int)- Maximum distance to left region.
        
        Returns:
            output - A list of named tuples containing five elements:
                chrom (str)- Name of chromosome.
                start (int)- Start of interval.
                end (int)- End of interval.
                strand (str)- Interval strand.
                data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        if strand == '*':
            strand = '+-.'
        self.check_position(chrom, position, strand)
        if not isinstance(dist, int):
            raise TypeError('dist not integer')
        if dist < 0:
            raise ValueError('dist negative')
        dist += 1
        #  Loop through strands and find nearest upstream interval
        minDist = None
        for s in strand:
            # Find left interval end and distance
            before = self.intervals[(chrom, s)].before(
                position + 1, max_dist = dist)
            try:
                end = before[0].end
            except IndexError:
                continue
            # Store minimum distance and strand
            distance = position - end
            if minDist is None or distance < minDist:
                minDist = distance
                minStrand = s
                minEnd = end
            elif distance == minDist:
                minStrand += s
        # Find upstream intervals
        output = []
        if not minDist is None:
            for s in minStrand:
                output.extend(self.find_end(chrom, minEnd, s))
            if len(minStrand) > 1:
                output = sorted(output, key=lambda x: x[1:4])
        # Return data
        return(output, minDist)
    
    def find_right(
            self, chrom, position, strand = '*', dist=1000000000
        ):
        ''' A function to return data with leftmost intervals
        
        Args:
            chrom (str)- Name of chromosome.
            position (int)- Position on chromosome.
            strand (str)- Must be a combination of '+/-/.' or '*'.
            dist (int)- Maximum distance to right region.
        
        Returns:
            output - A list of named tuples containing five elements:
                chrom (str)- Name of chromosome.
                start (int)- Start of interval.
                end (int)- End of interval.
                strand (str)- Interval strand.
                data (dict)- Dictionary of interval decriptors.
        
        '''
        # Check arguments
        if strand == '*':
            strand = '+-.'
        self.check_position(chrom, position, strand)
        if not isinstance(dist, int):
            raise TypeError('dist not integer')
        if dist < 0:
            raise ValueError('dist negative')
        dist += 1
        #  Loop through strands and find nearest upstream interval
        minDist = None
        for s in strand:
            # Find left interval end and distance
            after = self.intervals[(chrom, s)].after(
                position - 1, max_dist = dist)
            try:
                start = after[0].start
            except IndexError:
                continue
            # Store minimum distance and strand
            distance = start - position
            if minDist is None or distance < minDist:
                minDist = distance
                minStrand = s
                minStart = start
            elif distance == minDist:
                minStrand += s
        # Find upstream intervals
        output = []
        if not minDist is None:
            for s in minStrand:
                output.extend(self.find_start(chrom, minStart, s))
            if len(minStrand) > 1:
                output = sorted(output, key=lambda x: x[1:4])
        # Return data
        return(output, minDist)

# Perform unittesting
if __name__ == "__main__":
    
    import unittest
    
    class FindIntervalTest(unittest.TestCase):
        
        def setUp(self):
            # Create intervals for unittest
            inttup = collections.namedtuple(
                'inttup', ['chrom', 'start', 'end', 'strand', 'data'])
            self.int1 = inttup('chr1', 30, 49, '.', {'name':'A'})
            self.int2 = inttup('chr1', 30, 50, '+', {'name':'B'})
            self.int3 = inttup('chr1', 35, 55, '-', {'name':'C'})
            self.int4 = inttup('chr1', 40, 60, '.', {'name':'C'})
            self.int5 = inttup('chr1', 45, 65, '+', {'name':'D'})
            self.int6 = inttup('chr1', 50, 70, '-', {'name':'E'})
            self.int7 = inttup('chr1', 55, 75, '.', {'name':'F'})
            self.int8 = inttup('chr1', 60, 80, '+', {'name':'G'})
            self.int9 = inttup('chr1', 65, 85, '-', {'name':'H'})
            self.int10 = inttup('chr1', 70, 90, '.', {'name':'I'})
            self.int11 = inttup('chr1', 71, 90, '+', {'name':'J'})
            # Build interval object for testing
            self.fi = FindInterval(chroms={'chr1':120})
            for interval in (
                    self.int1, self.int2, self.int3, self.int4, self.int5,
                    self.int6, self.int7, self.int8, self.int9, self.int10,
                    self.int11
                ):
                self.fi.add_interval(*interval)
            # Function to compare list of tuples
            def tup_equal(list1, list2):
                if len(list1) != len(list2):
                    return(False)
                for tup1, tup2 in zip(list1, list2):
                    if tup1.chrom != tup2.chrom:
                        return(False)
                    if tup1.start != tup2.start:
                        return(False)
                    if tup1.end != tup2.end:
                        return(False)
                    if tup1.strand != tup2.strand:
                        return(False)
                    if tup1.data != tup2.data:
                        return(False)
                return(True)
            self.tup_equal = tup_equal

    class FindRightTest(FindIntervalTest):
        ''' Test FindInterval.find_left function '''
        
        def test_multistrand_right_1(self):
            right, dist = self.fi.find_right('chr1', 30, '*')
            expectedRight, expectedDist = [self.int1, self.int2], 0
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_multistrand_right_2(self):
            right, dist = self.fi.find_right('chr1', 15, '+.', dist=15)
            expectedRight, expectedDist = [self.int1, self.int2], 15
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_multistrand_right_3(self):
            right, dist = self.fi.find_right('chr1', 15, '+-', dist=14)
            expectedRight, expectedDist = [], None
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_right_1(self):
            right, dist = self.fi.find_right('chr1', 36, '-', dist=14)
            expectedRight, expectedDist = [self.int6], 14
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_right_2(self):
            right, dist = self.fi.find_right('chr1', 35, '-', dist=15)
            expectedRight, expectedDist = [self.int3], 0
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_right_3(self):
            right, dist = self.fi.find_right('chr1', 66, '-')
            expectedRight, expectedDist = [], None
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_right_1(self):
            right, dist = self.fi.find_right('chr1', 30, '+', dist=0)
            expectedRight, expectedDist = [self.int2], 0
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_right_2(self):
            right, dist = self.fi.find_right('chr1', 31, '+', dist=14)
            expectedRight, expectedDist = [self.int5], 14
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_right_3(self):
            right, dist = self.fi.find_left('chr1', 31, '+', dist=13)
            expectedRight, expectedDist = [], None
            self.assertTrue(self.tup_equal(right, expectedRight))
            self.assertTrue(dist == expectedDist)
    
    class FindLeftTest(FindIntervalTest):
        ''' Test FindInterval.find_left function '''
        
        def test_multistrand_left_1(self):
            left, dist = self.fi.find_left('chr1', 90)
            expectedLeft, expectedDist = [self.int10, self.int11], 0
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_multistrand_left_2(self):
            left, dist = self.fi.find_left('chr1', 90, '+.')
            expectedLeft, expectedDist = [self.int10, self.int11], 0
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_multistrand_left_3(self):
            left, dist = self.fi.find_left('chr1', 100, '+-')
            expectedLeft, expectedDist = [self.int11], 10
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_left_1(self):
            left, dist = self.fi.find_left('chr1', 100, '-', dist=15)
            expectedLeft, expectedDist = [self.int9], 15
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_left_2(self):
            left, dist = self.fi.find_left('chr1', 100, '-', dist=14)
            expectedLeft, expectedDist = [], None
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_negstrand_left_3(self):
            left, dist = self.fi.find_left('chr1', 55, '-')
            expectedLeft, expectedDist = [self.int3], 0
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_left_1(self):
            left, dist = self.fi.find_left('chr1', 65, '+', dist=0)
            expectedLeft, expectedDist = [self.int5], 0
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_left_2(self):
            left, dist = self.fi.find_left('chr1', 64, '+', dist=14)
            expectedLeft, expectedDist = [self.int2], 14
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
        
        def test_plusstrand_left_3(self):
            left, dist = self.fi.find_left('chr1', 64, '+', dist=13)
            expectedLeft, expectedDist = [], None
            self.assertTrue(self.tup_equal(left, expectedLeft))
            self.assertTrue(dist == expectedDist)
    
    class FindStartTest(FindIntervalTest):
        ''' Test FindInterval.find_start function '''
        
        def test_multistranded_start_1(self):
            start = self.fi.find_start('chr1', 30, '*')
            expected = [self.int1, self.int2]
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_multistranded_start_2(self):
            start = self.fi.find_start('chr1', 30, '+.')
            expected = [self.int1, self.int2]
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_multistranded_start_3(self):
            start = self.fi.find_start('chr1', 29, '+-.')
            expected = []
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_unstranded_start_1(self):
            start = self.fi.find_start('chr1', 40, '.')
            expected = [self.int4]
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_unstranded_start_2(self):
            start = self.fi.find_start('chr1', 41, '.')
            expected = []
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_unstranded_start_3(self):
            start = self.fi.find_start('chr1', 45, '.')
            expected = []
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_plusstrand_start_1(self):
            start = self.fi.find_start('chr1', 45, '+')
            expected = [self.int5]
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_plusstrand_start_2(self):
            start = self.fi.find_start('chr1', 44, '+')
            expected = []
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_negstrand_start_1(self):
            start = self.fi.find_start('chr1', 50, '-')
            expected = [self.int6]
            self.assertTrue(self.tup_equal(start, expected))
        
        def test_negstrand_start_2(self):
            start = self.fi.find_start('chr1', 51, '-')
            expected = []
            self.assertTrue(self.tup_equal(start, expected))

    class FindEndTest(FindIntervalTest):
        ''' Test FindInterval.find_end function '''
        
        def test_multistranded_end_1(self):
            end = self.fi.find_end('chr1', 90, '*')
            expected = [self.int10, self.int11]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_multistranded_end_2(self):
            end = self.fi.find_end('chr1', 90, '.+')
            expected = [self.int10, self.int11]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_multistranded_end_3(self):
            end = self.fi.find_end('chr1', 75, '+-')
            expected = []
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_unstranded_end_1(self):
            end = self.fi.find_end('chr1', 90, '.')
            expected = [self.int10]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_unstranded_end_2(self):
            end = self.fi.find_end('chr1', 60, '.')
            expected = [self.int4]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_unstranded_end_3(self):
            end = self.fi.find_end('chr1', 70, '.')
            expected = []
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_plusstrand_end_1(self):
            end = self.fi.find_end('chr1', 90, '+')
            expected = [self.int11]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_plusstrand_end_2(self):
            end = self.fi.find_end('chr1', 65, '+')
            expected = [self.int5]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_negstrand_end_1(self):
            end = self.fi.find_end('chr1', 85, '-')
            expected = [self.int9]
            self.assertTrue(self.tup_equal(end, expected))
        
        def test_negstrand_end_2(self):
            end = self.fi.find_end('chr1', 84, '-')
            expected = []
            self.assertTrue(self.tup_equal(end, expected))

    class FindOverlapTest(FindIntervalTest):
        ''' Test FindInterval.find_overlap function '''
            
        def test_allstrand_overlap_1(self):
            overlap = self.fi.find_overlap('chr1', 50, 70, '*')
            expected = [self.int3, self.int4, self.int5, self.int6, self.int7,
                self.int8, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_allstrand_overlap_2(self):
            overlap = self.fi.find_overlap('chr1', 49, 71, '*')
            expected = [self.int2, self.int3, self.int4, self.int5, self.int6,
                self.int7, self.int8, self.int9, self.int10]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_allstrand_overlap_3(self):
            overlap = self.fi.find_overlap('chr1', 0, 30, '*')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_plusstrand_overlap_1(self):
            overlap = self.fi.find_overlap('chr1', 50, 71, '+')
            expected = [self.int5, self.int8]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_plustrand_overlap_2(self):
            overlap = self.fi.find_overlap('chr1', 49, 72, '+')
            expected = [self.int2, self.int5, self.int8, self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_plustrand_overlap_3(self):
            overlap = self.fi.find_overlap('chr1', 0, 30, '+')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_negstrand_overlap_1(self):
            overlap = self.fi.find_overlap('chr1', 55, 65, '-')
            expected = [self.int6]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_negstrand_overlap_2(self):
            overlap = self.fi.find_overlap('chr1', 54, 66, '-')
            expected = [self.int3, self.int6, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_negstrand_overlap_3(self):
            overlap = self.fi.find_overlap('chr1', 0, 35, '-')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_unknownstrand_overlap_1(self):
            overlap = self.fi.find_overlap('chr1', 49, 70, '.')
            expected = [self.int4, self.int7]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_unknownstrand_overlap_2(self):
            overlap = self.fi.find_overlap('chr1', 48, 71, '.')
            expected = [self.int1, self.int4, self.int7, self.int10]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_unknownstrand_overlap_3(self):
            overlap = self.fi.find_overlap('chr1', 0, 30, '.')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_multistrand_overlap_1(self):
            overlap = self.fi.find_overlap('chr1', 50, 71, '-+')
            expected = [self.int3, self.int5, self.int6, self.int8, self.int9]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_multistrand_overlap_2(self):
            overlap = self.fi.find_overlap('chr1', 49, 72, '-+')
            expected = [self.int2, self.int3, self.int5, self.int6, self.int8,
                self.int9, self.int11]
            self.assertTrue(self.tup_equal(overlap, expected))
        
        def test_multistrand_overlap_3(self):
            overlap = self.fi.find_overlap('chr1', 90, 120, '+-')
            expected = []
            self.assertTrue(self.tup_equal(overlap, expected))
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindOverlapTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindStartTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindEndTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindLeftTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(FindRightTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
