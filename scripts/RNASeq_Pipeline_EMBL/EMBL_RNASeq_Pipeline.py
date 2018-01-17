#!/bin/python2.7
import argparse
import collections
import os
import re
import sys
scriptdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scriptdir)
import EMBL_RNASeq_Functions as rsf
import EMBL_RNASeq_Slurm as slurm

# Extract command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    'samples', type=str, help='file containing sample data')
parser.add_argument(
    'parameters', type=str, help='file containing parameter data')
parser.add_argument(
    'modules', type=str, help='file containing path data')
parser.add_argument(
    '--resubmit', help='resubmit failed samples', action='store_true')
args = parser.parse_args()
# Parse input files
samples = rsf.parseSampleFile(args.samples)
params = rsf.parseParameterFile(args.parameters)
paths, modules = rsf.parseModuleFile(args.modules)
# Create output directory
params['outdir'] = os.path.abspath(params['outdir'])
if not os.path.isdir(params['outdir']):
    os.mkdir(params['outdir'])
# Loop through samples and create slurm objects for each sample
slurmDict = collections.OrderedDict()
jobCounter = collections.Counter()
for sample, (prefix, indirs) in samples.items():
    # Create output file names and process existing files
    jobCounter['samples'] += 1
    outfiles = rsf.createOutFiles(params['outdir'], sample)
    sampleJobs = slurm.submitJobs()
    # Process prexisting jobs
    if os.path.isfile(outfiles['slurm']):
        # Extract jobstatus
        jobStatus = slurm.parseSlurmFile(outfiles['slurm'])
        # Clear outout directories containing failed jobs if resubmit requested
        if 'failed' in [x[2] for x in jobStatus] and args.resubmit:
            for d in (outfiles['outdir'], outfiles['seqcdir']):
                for f in os.listdir(d):
                    p = os.path.join(d, f)
                    if os.path.isfile(p):
                        os.remove(p)
        # Else simply print status of jobs and skip resubmission
        else:
            print(sample)
            for jobid, jobname, status in jobStatus:
                print('  {}\t{}\t{}'.format(jobname, jobid, status))
                jobCounter['jobs'] += 1
                jobCounter[status] += 1
            continue
    # Extract FASTQ files and concatenate if required
    read1, read2 = rsf.findFastq(prefix, indirs)
    # Concatenate reads and rename input FASTQ files
    catJobList = []
    if len(read1) > 1:
        cat1Command = 'zcat {} | gzip > {}'.format(
            ' '.join(read1), outfiles['cat1'])
        cat1JobID = sampleJobs.add(cat1Command, name=sample + '.zcat1')
        catJobList.append(cat1JobID)
    elif len(read1) == 1:
        if os.path.islink(outfiles['cat1']):
            os.remove(outfiles['cat1'])
        os.symlink(read1[0], outfiles['cat1'])
    if len(read2) > 1:
        cat2Command = 'zcat {} | gzip > {}'.format(
            ' '.join(read2), outfiles['cat2'])
        cat2JobID = sampleJobs.add(cat2Command, name=sample + '.zcat2')
        catJobList.append(cat2JobID)
    elif len(read2) == 1:
        if os.path.islink(outfiles['cat2']):
            os.remove(outfiles['cat2'])
        os.symlink(read2[0], outfiles['cat2'])
    else:
        outfiles['cat2'] = None
        outfiles['trim2'] = None
    # Set single end variable
    if outfiles['trim2'] is None:
        singleEnd = True
    else:
        singleEnd = False
    # Perform read trimming
    cutadaptCommand = rsf.cutadapt(
        read1In=outfiles['cat1'], read1Out=outfiles['trim1'],
        read2In=outfiles['cat2'], read2Out=outfiles['trim2'],
        quality=params['basequal'], adapter=params['adapter'],
        length=params['minlength'], path=paths['cutadapt'],
        overlap=params['overlap'], error=params['error']
    )
    cutadaptJobID = sampleJobs.add(
        cutadaptCommand, processors=1, memory=6, depend=catJobList,
        modules=modules['cutadapt'], stdout=outfiles['trimlog'],
        stderr=outfiles['trimlog'], name=sample + '.cutadapt'
    )
    # Perform fastqc
    fastqcCommand = []
    for fastq in (
            outfiles['cat1'], outfiles['cat2'], outfiles['trim1'],
            outfiles['trim2']
        ):
        if fastq is None:
            continue
        fastqcCommand.append(
            rsf.fastQC(
                inFile=fastq, outDir=outfiles['outdir'],
                path=paths['fastqc']))
    fastqcCommand = '; '.join(fastqcCommand)
    # Submit FASTQC commands
    fastqcJobID = sampleJobs.add(
        fastqcCommand, processors=1, memory=4, depend=[cutadaptJobID],
        modules=modules['fastqc'], stderr=outfiles['fastqclog'],
        stdout=outfiles['fastqclog'], name=sample + '.fastqc'
    )
    # Perform alignment with STAR
    starCommand = rsf.starAlign(
        indexDir=params['index'], outPrefix=outfiles['prefix'],
        read1=outfiles['trim1'], read2=outfiles['trim2'],
        path=paths['star'], threads=16, rg='1', lb=prefix, sm=sample,
        pl='illumina'
    )
    starJobID = sampleJobs.add(
        starCommand, processors=16, memory=6, depend=[cutadaptJobID],
        modules=modules['star'], name=sample + '.star',
        stdout=outfiles['starlog'], stderr=outfiles['starlog']
    )
    # Perform bam sort and idexing
    sortCommand = rsf.bamsort(
        inFile=outfiles['starbam'], outFile=outfiles['sortbam'], threads=8,
        memory=4, path=paths['samtools']
    )
    sortJobID = sampleJobs.add(
        sortCommand, processors=8, memory=6, depend=[starJobID],
        modules=modules['samtools'], name=sample + '.sort',
        stdout=outfiles['sortlog'], stderr=outfiles['sortlog']
    )
    # Mark duplicates in BAM
    mdupCommand = rsf.markDuplicates(
        inBam=outfiles['sortbam'], outBam=outfiles['mdupbam'],
        logFile=outfiles['mduplog1'], picardPath=paths['picard'],
        memory=30
    )
    mdupJobID = sampleJobs.add(
        mdupCommand, depend=[sortJobID], modules=['picard'],
        stdout=outfiles['mduplog2'], stderr=outfiles['mduplog2'], processors=6,
        memory=6, name=sample + '.markdup'
    )
    # Extract htseq gene counts
    htseqCommand = rsf.htseq(
        bam=outfiles['mdupbam'], gtf=params['gtf'], path=paths['htseq'],
        feature=params['feature'], attrid=params['id'], mode=params['mode'],
        stranded=params['stranded'], mapq=params['mapq']
    )
    htseqJobID = sampleJobs.add(
        htseqCommand, depend=[mdupJobID], modules=modules['htseq'],
        stdout=outfiles['genecounts'], stderr=outfiles['htseqlog'],
        processors=1, memory=6, name=sample + '.htseq'
    )
    # Collect picard RNA-Seq metrics
    metricCommand = rsf.rnaseqMetric(
        bam=outfiles['mdupbam'], output=outfiles['metrlog1'],
        refflat=params['refflat'], strand=params['strand'],
        rrna=params['rrna'], path=paths['picard'], memory=30
    )
    metricJobID = sampleJobs.add(
        metricCommand, stdout=outfiles['metrlog2'], name=sample + '.rnametric',
        stderr=outfiles['metrlog2'], depend=[mdupJobID],
        modules=modules['picard'], memory=6, processors=6
    )
    # Collect picard alignment metrics
    alsumCommand = rsf.alignMetrics(
        bam=outfiles['mdupbam'], output=outfiles['alsumlog1'],
        fasta=params['fasta'], path=paths['picard'], memory=30
    )
    alsumJobID = sampleJobs.add(
        alsumCommand, stdout=outfiles['alsumlog2'],
        stderr=outfiles['alsumlog2'], modules=modules['picard'], memory=6,
        processors=6, depend=[mdupJobID], name=sample + '.alignsum'
    )
    # Store commands
    slurmDict[sample] = sampleJobs

###############################################################################
## Submit jobs and print metrics
###############################################################################
for sample, sampleJobs in slurmDict.items():
    # Open output file
    slurmfile = os.path.join(params['outdir'], sample, sample + '.slurm')
    with open(slurmfile, 'w') as sf:
        # Submit jobs and write output to file
        jobList = sampleJobs.submit(check_sub=False)
        outstring = sampleJobs.output(jobList)
        sf.write(outstring)
    # Write data to stdout
    print(sample)
    for job in jobList:
        jobid, jobname = job[:2]
        print('  {}\t{}\tsubmitted'.format(jobname, jobid))
        jobCounter['jobs'] += 1
        jobCounter['submitted'] += 1
# Print output numbers
print('All Job Summary')
print('  {:<12}{}'.format('samples', jobCounter['samples']))
print('  {:<12}{}'.format('jobs', jobCounter['jobs']))
print('  {:<12}{}'.format('submitted', jobCounter['submitted']))
print('  {:<12}{}'.format('pending', jobCounter['pending']))
print('  {:<12}{}'.format('running', jobCounter['running']))
print('  {:<12}{}'.format('failed ', jobCounter['failed']))
print('  {:<12}{}'.format('completed', jobCounter['completed']))
