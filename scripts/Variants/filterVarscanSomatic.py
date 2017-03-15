'''filterVarscanSomatic.py
    
    Usage:
    
    varscan.py <varscandata> <controlbam> <testbam> <outfile> <parameters>
    
'''
# Import modules
import collections
import os
from ngs_python.bam.pysam_variants import PysamVariants
from general_python import docopt
args = docopt.docopt(__doc__,version = 'v1')
# Parse parameter file
params = collections.OrderedDict()
with open(args['<parameters>']) as infile:
    for line in infile:
        linelist = line.strip().split('\t')
        if linelist[0] in params:
            print(linelist)
            print(params)
            raise ValueError('repeat parameter: {}'.format(linelist[0]))
        params[linelist[0]] = linelist[1]
# Create function to return boolean for parameters
def return_bool(arg):
    if arg.startswith('t') or arg.startswith('T'):
        return(True) 
    elif arg.startswith('f') or arg.startswith('F'):
        return(False)
    else:
        raise ValueError('could not convert {} to bool'.format(arg))
# Check parameters
params['min base quality'] = int(params['min base quality'])
if not 0 < params['min base quality'] <= 255:
    raise ValueError('min base quality must be >0 and <=256')
params['min map quality'] = int(params['min map quality'])
if not 0 <= params['min map quality'] <= 255:
    raise ValueError('min map quality must be >=0 and <=255')
params['mean map quality'] = int(params['mean map quality'])
if not 0 <= params['mean map quality'] <= 255:
    raise ValueError('mean map quality must be >=0 and <=255')
params['min control coverage'] = int(params['min control coverage'])
if params['min control coverage'] < 1:
    raise ValueError('min control coverage must be >0')
params['max control frequency'] = float(params['max control frequency'])
if not 0 <= params['max control frequency'] <= 1:
    raise ValueError('max control frequency must be >0 and <=1')
params['min test coverage'] = int(params['min test coverage'])
if params['min test coverage'] < 1:
    raise ValueError('min test coverage must be >0')
params['min test frequency'] = float(params['min test frequency'])
if not 0 < params['min test frequency'] <= 1:
    raise ValueError('min test frequency must be >0 and <=1')
params['keep duplicates'] = return_bool(params['keep duplicates'])
params['keep secondary'] = return_bool(params['keep secondary'])
params['keep supplementary'] = return_bool(params['keep supplementary'])
params['keep qc fail'] = return_bool(params['keep qc fail'])
params['proper pairs'] = return_bool(params['proper pairs'])
# Read in varscan file
varlist = []
with open(args['<varscandata>']) as infile:
    header = infile.next().strip().split('\t')
    stateIndex = header.index('somatic_status')
    chromIndex = header.index('chrom')
    posIndex = header.index('position')
    refIndex = header.index('ref')
    varIndex = header.index('var')
    for line in infile:
        linelist = line.strip().split('\t')
        if linelist[stateIndex] == 'Somatic':
            chrom = linelist[chromIndex]
            pos = int(linelist[posIndex])
            ref = linelist[refIndex]
            var = linelist[varIndex]
            varlist.append((chrom, pos, ref, var))
# Parse varscan file
with PysamVariants([args['<controlbam>'], args['<testbam>']]) as pv:
    variantData = pv.collect_variant_counts(
        variants=varlist, mapq=params['min map quality'],
        baseq=params['min base quality'], unmapped=False,
        qcfail=params['keep qc fail'], duplicate=params['keep duplicates'],
        secondary=params['keep secondary'],
        supplementary=params['keep supplementary'],
        properpair=params['proper pairs'], threads=2)
# Filter variant data
with open(args['<outfile>'], 'w') as outfile:
    # Print parameters
    outfile.write('# norm: {}\n'.format(os.path.abspath(args['<controlbam>'])))
    outfile.write('# test: {}\n'.format(os.path.abspath(args['<testbam>'])))
    for p in params:
        outfile.write('# {}: {}\n'.format(p, params[p]))
    # Filter and print variants
    outfile.write('chr\tpos\tref\tvar\tnorm_cov\tnorm_freq\tnorm_mapq'\
        '\ttest_cov\ttest_freq\tmapq\n')
    for count, variant in enumerate(varlist):
        # Extract and check control data
        control = variantData[0][count]
        controlCov = control[0] + control[2]
        if controlCov < params['min control coverage']:
            continue
        controlFreq = control[2] / float(controlCov)
        if controlFreq > params['max control frequency']:
            continue
        # Extract and check test data
        test = variantData[1][count]
        testCov = test[0] + test[2]
        if testCov < params['min test coverage']:
            continue
        testFreq = test[2] / float(testCov)
        if testFreq < params['min test frequency']:
            continue
        # Extract and check mapq data
        if min(control[4], test[4]) < params['mean map quality']:
            continue
        # Write acceptable data to file
        outfile.write('{}\t{}\t{:.3f}\t{}\t{}\t{:.3f}\t{}\n'.format(
            '\t'.join(map(str, variant)), controlCov, controlFreq, control[4],
            testCov, testFreq, test[4]))
