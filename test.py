import multiprocessing
import sys
import os
sys.path.append('/home/adam/github/')
from ngs_analysis.bam import bamExtract
print 'Parent: ' + str(os.getpid())
# Create process to read in read1 FASTQ file
pipes1 = multiprocessing.Pipe(False)
process1 = multiprocessing.Process(
    target = bamExtract.readToPipe,
    args = ('read1_test.fastq.gz', pipes1)
)
process1.start()
print 'Read1: ' + str(process1.pid)
pipes1[1].close()
# Create process to read in read2 FASTQ file
pipes2 = multiprocessing.Pipe(False)
process2 = multiprocessing.Process(
    target = bamExtract.readToPipe,
    args = ('read2_test.fastq.gz', pipes2)
)
process2.start()
print 'Read2: ' + str(process2.pid)
pipes2[1].close()
# Create process to write to output FASTQ file
pipes3 = multiprocessing.Pipe(False)
process3 = multiprocessing.Process(
    target = bamExtract.writeFromPipe,
    args = ('readOut_test.fastq.gz', pipes3)
)
process3.start()
print 'Write: '+ str(process3.pid)
count = 0
while True:
    try:
        read1 = pipes1[0].recv()
        read2 = pipes2[0].recv()
    except EOFError:
        print read1
        print read2
        break
    count += 2
    if not (count % 100000):
        print str(count) + ' read'
    pipes3[1].send(read1)
    pipes3[1].send(read2)
pipes3[1].close()
process1.join()
print 'joined 1'
process2.join()
print 'joined 2'
process3.join()
print 'joined 3'
