import collections

class bedData(object):
    
    def __init__(self, chromList = None):
        self.bedDict = collection.OrderedDict():
        if chromList:
            for chrom in self.bedDict[chrom] = ''
    
    def add(self, *args):
        # Unpack arguments
        if len(args) == 3:
            chrom, start, end = args
        elif len(args) == 1:
            chrom, start, end = args[0]
        else:
            raise ValueError('3 args, or a tuple/list of 3 args required')
        # Check chromosome is a string
        if not isinstance(chrom, str):
            print 'Invalid chromosome', chrom
            raise ValueError('First arg must be a chromosome string')
        # Check start is integer
        if not isinstance(start, int):
            print 'Invalid start', start
            raise ValueError('Second arg must be a start integer')
        # Check end is integer
        if not isinstance(end, int):
            print 'Invalid end', end
            raise ValueError('Third arg must be an end integer')
        # Check interval
        if not start >= 0 and end > start:
            print 'Invalid interval', start, end
            raise ValueError('Start must be >= 0 and less than end')
        # Add data to dictionary
        if chrom in self.bedDict:
            self.bedDict[chrom].append((start, end))
        else:
            self.bedDict[chrom] = [(start, end)]
    
    def sortStart(self):
        # Loop through chromsomes
        for chrom in self.bedDict:
            # Sort intervals by start
            self.bedDict[chrom].sort()
    
    def mergeIntervals(self):
        # Create variable to store output
        mergeDict = collections.OrderedDict()
        # Sort intervals
        self.sortStart()
        # Loop through chrosome
        for chrom in self.bedDict:
            # Create variables
            result = []
            currentStart = -1
            currentStop = -1
            # Loop through sorted intervals
            for start, stop in self.bedDict[chrom]:
                if start > currentStop:
                    result.append((start, stop))
                    currentStart, currentStop = start, stop
                else:
                    result[-1] = (currentStart, stop)
                    currentStop = max(currentStop, stop)
            # Store merged intervals
            self.mergeDict[chrom] = result
        # Return merged intervals
        return(mergeDict)
    
    def bedGraph(self):
        # Create variable to store output
        bedGraph = collections.OrderedDict()
        # Loop through chromsomes
        for chrom in self.bedDict:
            # Create variables
            countDict = collections.defaultdict(int)
            # Loop through intervals and add counts to dictionary
            for start, stop in self.bedDict[chrom]:
                # Add counts
                for position in xrange(start, stop):
                    countDict[position] += 1
            # Sequenctially process genomic positions
            for count, position in enumerate(sorted(countDict.keys())):
                # Extract count
                coverage = countDict[position]
                # Add first interval
                if not count:
                    intervalList = [[position, position, coverage]]
                # Process concurrent intervals
                if (position == intervalList[-1][1] and
                    coverage == intervalList[-1][2]):
                    intervalList[-1][1] += 1
                # Else store interval
                else:
                    intervalList.append([position, position + 1,
                        coverage])
            # Add final interval
            bedGraph[chrom] = intervalList
        # Return data
        return(bedGraph)

bd = bedData()
bd.add(('chr1',20,50))
bd.add('chr1',10,20)
print bd.bedDict
print bd.bedGraph()
