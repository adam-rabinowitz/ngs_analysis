import pysam
import collections
import multiprocessing
import functools
import subprocess
import pandas as pd
from ngs_python.bam.pysam_bam import PysamBAM

class PysamVariants(PysamBAM):
    
    def _variant_basecalls(
            self, bam, filterfunc, chrom, pos, ref, var
        ):
        # Create base tuple for filtering base calls
        outlist = []
        basetuple = (ref, var)
        strandDict = {False : '+', True : '-'}
        # Loop through pileup column and extract basecalls
        for read, index in self._getpileup(
                bam=bam, filterfunc=filterfunc, chrom=chrom, pos=pos,
                mapq=None
            ):
            # Extract and remove unwanted bases
            base = read.query_sequence[index]
            if base not in basetuple:
                continue
            baseq = read.query_qualities[index]
            mapq = read.mapping_quality
            strand = strandDict[read.is_reverse]
            outlist.append((base, baseq, mapq, strand))
        # Return data
        return(outlist)
    
    def _variant_counts(
            self, bamindex, filterfunc, variants, mapq, baseq
        ):
        # Create output variable and extract variant data from indexed bam
        outlist = []
        with pysam.AlignmentFile(self.bampaths[bamindex], 'rb') as bam:
            for chrom, pos, ref, var in variants:
                basecalls = self._variant_basecalls(
                    bam=bam, filterfunc=filterfunc, chrom=chrom, pos=pos,
                    ref=ref, var=var)
                # Extract base calls for filtered reads
                count = collections.defaultdict(int)
                for read in basecalls:
                    # Skip low mapping quality and base qualities
                    count['mapq_sum'] += read[2]
                    count['mapq_len'] += 1
                    if read[2] < mapq:
                        continue
                    if read[1] < baseq:
                        continue
                    # Extract reference and variant counts
                    if read[0] == ref:
                        count['ref'] += 1
                        if read[3] == '+':
                            count['ref_for'] += 1
                    elif read[0] == var:
                        count['var'] += 1
                        if read[3] == '+':
                            count['var_for'] += 1
                # Add variant data to output
                meanmapq = count['mapq_sum'] / count['mapq_len']
                outtuple = (
                    count['ref'], count['ref_for'], count['var'],
                    count['var_for'], meanmapq)
                outlist.append(outtuple)
        # Return data
        return(outlist)
        
    def _variant_counts_process(
            self, inqueue, outqueue, filterfunc, variants, mapq, baseq
        ):
        # Extract indices from queue
        while True:
            # Extract and check index
            index = inqueue.get()
            if index is None:
                inqueue.close()
                outqueue.close()
                break
            # Extract data and add to counts
            variantdata = self._variant_counts(
                bamindex=index, filterfunc=filterfunc, variants=variants,
                mapq=mapq, baseq=baseq)
            outqueue.put((index, variantdata))
    
    def collect_variant_counts(
            self, variants, mapq = 20, baseq = 20, unmapped = False,
            qcfail = False, duplicate = False, secondary = False,
            supplementary = False, properpair=True, threads = None
        ):
        ''' Function to extract variant metrics for bam file in object.
        
        Args:
            variants - Iterable returning chromosome, position, reference base
                and variant base of variants.
            mapq (int)- Minimum alignment mapping quality.
            baseq (int)- Minimum alignment base quality.
            umapped (bool)- Include unmapped reads.
            qcfail (bool)- Include qc failed reads.
            duplicate (bool)- Include duplicate alignments.
            secondary (bool)- Include secondary alignments.
            supplementary (bool)- Include supplementary alignments.
            propepair (bool)- Require properly paired alignments.
            threads (int)- Number of threads.
        
        Returns:
            outlist - A list, of lists of tuples. The top level list has a list
                of tuples, one for each bam file. The second level is list of
                tuples, one for each variant. The third level is a five
                element tuple consisiting of the following five elements:
                1) filtered count of reference bases.
                2) number of reference bases on the forward strand.
                3) filtered count of variant bases.
                4) number of variant bases on the forward strand.
                5) mean mapping quality of all reference and variant bases.
        
        '''
        # Process thread argument
        if threads is None:
            threads = len(self.bampaths)
        else:
            if not isinstance(threads, int):
                raise TypeError('threads must be integer')
            if not 1 <= threads <= 16:
                raise ValueError('threads must be >=1 and <= 16')
        # Create flag for filtering reads
        filterfunc = self._flagfilter(
            unmapped=unmapped, qcfail=qcfail, duplicate=duplicate,
            secondary=secondary, supplementary=supplementary,
            properpair=properpair)
        # Sequentially process BAM files
        if threads == 1 or len(self.bampaths) == 1:
            outlist = [None] * len(self.bampaths)
            for index in range(len(self.bampaths)):
                vardata = self._variant_counts(
                    bamindex=index, filterfunc=filterfunc, variants=variants,
                    mapq=mapq, baseq=baseq)
                outlist[index] = vardata
            return(outlist)
        # Parallel process BAM files
        else:
            # Create queues
            inqueue = multiprocessing.Queue()
            outqueue = multiprocessing.Queue()
            # Create processes
            processlist = []
            for x in range(threads):
                process = multiprocessing.Process(
                    target=self._variant_counts_process,
                    args=(inqueue, outqueue, filterfunc, variants, mapq,
                        baseq))
                process.start()
                processlist.append(process)
            # Add bam indices to inqueue and terminating None values
            for x in range(len(self.bampaths)):
                inqueue.put(x)
            for x in range(threads):
                inqueue.put(None)
            # Collect data from outqueue
            outlist = [None] * len(self.bampaths)
            while True:
                if not None in outlist:
                    break
                index, vardata = outqueue.get()
                outlist[index] = vardata
            inqueue.close()
            outqueue.close()
            # Close processes and return data
            for process in processlist:
                process.join()
            return(outlist)

class variants(object):
    
    def __init__(self, bamList):
        # Check bams contain identical reference sequence
        self.bams = bamList
        for count, bam in enumerate(bamList):
            inbam = pysam.AlignmentFile(bam)
            length = collections.OrderedDict(
                zip(inbam.references, inbam.lengths))
            if count:
                if length != self.length:
                    raise ValueError('Chromosomes in BAM files vary')
            else:
                self.length = length
            inbam.close()
    
    def __varCountSNP(
            self, bam, connection, variantList, minMapQ=20, minBaseQ=20,
            rmDup=True, rmSec=True, rmSup=True
        ):
        ''' Function to extract specified reference and variant counts at
        specied regions within a BAM file.
        
        Args:
            bam (str)- Path to BAM file.
            variantList (list)- A list of four element tuples containing:
                chromosome, position, reference base and variant base.
        
        Returns:
            countTuple - A tuple of the following four elements: reference
                count, forward strand reference count, variant count, forward
                strand variant count.
        
        '''
        # Process connections
        connIn, connOut = connection
        connIn.close()
        # Open BAM file and create output
        inbam = pysam.AlignmentFile(bam)
        outList = [None] * len(variantList)
        # Create flagFilter
        flagFilter = 516
        if rmDup:
            flagFilter += 1024
        if rmSec:
            flagFilter += 256
        if rmSup:
            flagFilter += 2048
        # Loop through variants
        for count, (chrom, position, ref, var) in enumerate(variantList):
            # Create output counters
            refCount, refForCount, varCount, varForCount = 0, 0, 0, 0
            counter = collections.defaultdict(int)
            # Find column in pileup
            position -= 1
            for column in inbam.pileup(
                    chrom, position - 1, position + 1, stepper='samtools'
                ):
                # Loop through reads in desired column
                if column.pos == position:
                    for read in column.pileups:
                        # Skip reads without base calls
                        if read.is_del or read.is_refskip:
                            continue
                        # Extract query data
                        alignment = read.alignment
                        query_position = read.query_position
                        # Skip gapped reads
                        if alignment.flag & flagFilter:
                            continue
                        # Check base and mapping quality
                        mapQ = alignment.mapping_quality
                        if mapQ < minMapQ:
                            continue
                        baseQ = alignment.query_qualities[query_position]
                        counter[baseQ] += 1
                        if baseQ < minBaseQ:
                            continue
                        # Extract base call
                        base = alignment.query_sequence[query_position]
                        if base == ref:
                            refCount += 1
                            if not alignment.is_reverse:
                                refForCount += 1
                        if base == var:
                            varCount += 1
                            if not alignment.is_reverse:
                                varForCount += 1
            # Add data to output
            outList[count] = (refCount, refForCount, varCount, varForCount)
        # Close input BAM and return data
        inbam.close()
        varnames = [':'.join(map(str,x)) for x in variantList]
        outDF = pd.DataFrame(
            outList,
            index = [':'.join(map(str,x)) for x in variantList],
            columns = ['Ref', 'RefFor', 'Var', 'VarFor']
        )
        connOut.send(outDF)
        connOut.close()
    
    def varCountSNP(
            self, variantList, minMapQ=20, minBaseQ=20, rmDup=True, rmSec=True,
            rmSup=True
        ):
        # Create partial function
        partialFunc = functools.partial(
            self.__varCountSNP, variantList=variantList, minMapQ=minMapQ,
            minBaseQ=minBaseQ, rmDup=rmDup, rmSec=rmSec,
            rmSup=rmSup
        )
        # Create an store processes
        processList = []
        for bam in self.bams:
            connection = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target=partialFunc,
                args=(bam, connection)
            )
            process.start()
            connection[1].close()
            processList.append((process, connection[0]))
        # Extract data from processes
        outList = []
        for process, conn in processList:
            outList.append(conn.recv())
            conn.close()
            process.join()
        return(outList)
