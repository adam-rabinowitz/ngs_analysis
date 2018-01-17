import collections
import gzip
import os
import re
import sys
from ngs_python.intervals import find_intervals

def gtf_generator(
        gtf, fixedwidth=True
    ):
    ''' A generator returning parsed lines from a gtf file
    
    Args:
        gtf (str)- Path to gtf/gff file
        fixedwidth (bool)- If True generator will raise error if line has
            less than 9 elements. If False such lines are skipped.
    
    Yields:
        outtup - A named tuple containing the following 9 elements:
            seqname (str)- Sequence name
            source (str)- Source of annotation
            feature (str)- Feature type
            start (int)- Start position
            end (int)- End position
            score (int, None)- Either integer or None.
            strand (str)- Strand.
            frame (str)- Either 0, 1, 2 or None.
            attributes (dict)- Keys are atrributes and values are values.
    
    '''
    # Check arguments
    if not isinstance(fixedwidth, bool):
        raise TypeError('fixedwidth must be bool')
    # Create named tuple for returning data
    entry = collections.namedtuple('entry', [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
        'frame', 'attributes'])
    # Set function to open file
    if gtf.endswith('.gz'):
        openfunc = gzip.open
    else:
        openfunc = open
    # Open file and loop through lines
    with openfunc(gtf) as infile:
        for line in infile:
            # Skip comments
            if line.startswith('#'):
                continue
            # Extract data from line
            linedata = line.strip().split('\t')
            try:
                seqname, source, feature = linedata[:3]
                start, end = map(int, linedata[3:5])
                score, strand, frame, attributes = linedata[5:9]
            except ValueError:
                if fixedwidth:
                    raise ValueError('too many values to upnack')
                else:
                    continue
            # Process attributes
            attributes = attributes.rstrip(';')
            attributes = attributes.replace('"', '')
            attributes = re.split(';\s*', attributes)
            attributeDict = {}
            for a in attributes:
                if a.startswith('#'):
                    continue
                a = re.split('[ =]', a)
                attributeDict[a[0]] = a[1]
            # Create tuple and yield
            outtup = entry(
                seqname, source, feature, start, end, score, strand, frame,
                attributeDict)
            yield(outtup)

def find_tss(
        gtf, feature='exon', identifier='gene_id', fixedwidth=True,
        skipanom=False
    ):
    ''' Function to parse gtf/gff files to extract position of the
    transcriptional start site of transcripts/genes based on the location
    of the defined features.
    
    Args:
        gtf (str)- Path to GTF/GFF file.
        feature (str)- GTF/GFF feature used to identify transcriptional
            start sites.
        identifier (str)- Attribute used to group features.
        fixedwidth (bool)- If True  will raise error if line in gtf/ggf
            has less than 9 elements. If False such lines are skipped.
        skipanom (bool)- Skip anomalous elements with ambiguos 
    
    '''
    # Create dictionary to store data and loop through gtf file
    featureDict = collections.defaultdict(list)
    tssDict = {}
    for entry in gtf_generator(
            gtf, fixedwidth=fixedwidth
        ):
        # Skip non exonic features
        if entry.feature != feature:
            continue
        # Extract feature id and store data
        featureID = entry.attributes[identifier]
        featureDict[featureID].append((
            entry.seqname, entry.start, entry.end, entry.strand))
    # Process features
    for featureID, featureData in featureDict.items():
        # Set intial values for location
        chrom = featureData[0][0]
        strand = featureData[0][3]
        if strand == '+':
            index = 1
            comparison = min
        elif strand == '-':
            index = 2
            comparison = max
        else:
            raise ValueError('Unexpected strand value')
        location = featureData[0][index]
        # Loop through and process additional data
        okay = True
        for fd in featureData[1:]:
            if fd[0] != chrom:
                error = '{} has ambiguous chromosome'.format(featureID)
                if skipanom:
                    sys.stderr.write(error + '\n')
                    okay = False
                    break
                else:
                    raise ValueError(error)
            if fd[3] != strand:
                error = '{} has ambiguous strand'.format(featureID)
                if skipanom:
                    sys.stderr.write(error + '\n')
                    okay = False
                    break
                else:
                    raise ValueError(error)
            location = comparison(location, fd[index])
        # Store data in output
        if okay:
            tssDict[featureID] = (chrom, location, strand)
    # Return data
    return(tssDict)

# Create transcript tuple class for storing transcript data
trantup = collections.namedtuple('trantup', ['tname', 'gname', 'gid'])

def extract_transcript_data(
        gtf, ttypes = ['rRNA', 'ncRNA', 'mRNA', 'snoRNA', 'pre_miRNA', 'snRNA',
        'tRNA', 'pseudogene'], fixedwidth=True
    ):
    ''' Function to extract transcript data for transcripts in flybase gff.
    Function parses all gene features and the following transcript features:
    rRNA, ncRNA, mRNA, snoRNA, pre-miRNA, snRNA, tRNA and pseudogenes.
    
    Args:
        gff (str)- Path to gff file.
        ttpes (list, tuple)- Lists acceptable transcript features
    
    Returns:
        tranDict (dict)- A Dictionary where the keys are the transcript IDs
            and the value is a named tuple containing the following elements:
            tname - transcript name
            gname - name of parent gene
            gid - id of parent gene
    
    '''
    # Check arguments
    if not os.path.isfile(gtf):
        raise IOError('could not find input file')
    if not isinstance(ttypes, (list, tuple)):
        raise TypeError('ttypes not a tuple or list')
    for t in ttypes:
        if not isinstance(t, str):
            raise TypeError('ttypes does contains non-strings')
    # Create dictionaries to store partial data
    id2nameGene = {}
    id2nameTran = {}
    child2parent = {}
    # Loop through gff and fill dictionaries
    for entry in gtf_generator(
            gtf, fixedwidth=fixedwidth, seqprefix=None, adjstart=0
        ):
        # Process gene feature
        if entry.feature == 'gene':
            # Extract gene features
            geneID = entry.attributes['ID']
            geneName = entry.attributes['Name']
            # Store gene features
            if geneID in id2nameGene:
                geneNamePrev = id2nameGene[geneID]
                if geneName != geneNamePrev:
                    error = 'gene id {} has two names: {} {}'.format(
                        geneID, geneNamePrev, geneName)
                    raise ValueError(error)
            else:
                id2nameGene[geneID] = geneName
        # Process transcript features
        if entry.feature in ttypes:
            # Extract transcript features
            tranID = entry.attributes['ID']
            tranName = entry.attributes['Name']
            geneID = entry.attributes['Parent']
            # Store transcript name
            if tranID in id2nameTran:
                tranNamePrev = id2nameGene[tranID]
                if tranName != tranNamePrev:
                    error = 'transcript id {} has two names: {} {}'.format(
                        tranID, tranNamePrev, tranName)
                    raise ValueError(error)
            else:
                id2nameTran[tranID] = tranName
            # Store transcript parent
            if tranID in child2parent:
                geneIDPrev = child2parent[tranID]
                if geneID != genIDPrev:
                    error = 'transcript {} has two genes: {} {}'.format(
                        tranID, geneIDPrev, geneID)
                    raise ValueError(error)
            else:
                child2parent[tranID] = geneID
    # Combine dictionaries
    tranDict = {}
    for tid in child2parent:
        tname = id2nameTran[tid]
        gid = child2parent[tid]
        gname = id2nameGene[gid]
        tidtup = trantup(tname, gname, gid)
        tranDict[tid] = tidtup
    for gid in id2nameGene:
        gname = id2nameGene[gid]
        gidtup = trantup(None, gname, gid)
        tranDict[gid] = gidtup
    return(tranDict)

def extract_feature_intervals(
        gff, features, chroms=None, seqprefix=None, adjstart=-1,
        fixedwidth=True
    ):
    ''' Function to build interval tree covering a specified gtf features
    
    Args:
        gtf (str)- Path to GTF/GFF file.
        features (dict)- A dictionary where the values are lists of
            features to extract from the input file and the key is the
            collective name of the features.
        chroms (dict)- A dictionary of chromosome lengths.
        adjstart (int)- Adjustment to make to interval start position.
            A value of -1 will convert 1-based closed intervals to 0-based
            semi-open intervals.
        seqprefix (str)- Prefix to add to sequence names.
    
    Returns:
        intervalDict (dict)- A dictionary where the values are FindInterval
            objects and the key is the collective name for the intervals.
    
    '''
    # Create variables for processing and storing output
    classDict = {}
    intervalDict = {}
    # Check feature dictionary
    if not isinstance(features, dict):
        raise TypeError('features not dictionary')
    for featureClass, featureList in features.items():
        if not isinstance(featureClass, str):
            raise TypeError('features key not string')
        if not isinstance(featureList, list):
            raise TypeError('features value not list')
        if not len(featureList) > 0:
            raise ValueError('features value an empty list')
    for featureClass, featureList in features.items():
        intervalDict[featureClass] = find_intervals.FindInterval()
        for feature in featureList:
            try:
                classDict[feature].append(featureClass)
            except KeyError:
                classDict[feature] = [featureClass]
    # Populate classDict and intervalDict
    for entry in gtf_generator(
            gff, seqprefix=seqprefix, adjstart=adjstart, fixedwidth=fixedwidth
        ):
        # Extract feature class
        try:
            featureClasses = classDict[entry.feature]
        except KeyError:
            continue
        # Store data
        for featureClass in featureClasses:
            intervalDict[featureClass].add_interval(
                chrom=entry.seqname, start=entry.start, end=entry.end,
                strand=entry.strand, data=entry.attributes)
    # Return data
    return(intervalDict)
:q
#def find_feature(
#        gtf, feature, chroms=None, adjstart=-1, seqprefix=None
#    ):
#    ''' Function to build interval tree covering a specified gtf feature
#    
#    Args:
#        gtf (str)- Path to GTF file.
#        feature (str)- Genome feature to find.
#        chroms (dict)- A dictionary of chromosome lengths.
#        adjstart (int)- Adjustment to make to interval start position.
#            A value of -1 will convert 1-based closed intervals to 0-based
#            semi-open intervals.
#        seqprefix (str)- Prefix to add to sequence names.
#    
#    '''
#    featuredata = find_intervals.FindInterval(chroms=chroms)
#    for entry in gtf_generator(
#            gtf, adjstart=adjstart, seqprefix=seqprefix
#        ):
#        if entry.feature != feature:
#            continue
#        featuredata.add_interval(
#            chrom=entry.seqname, start=entry.start, end=entry.end, strand=entry.strand,
#            data=entry.attributes)
#    return(featuredata)
#
#
#def find_features(
#        gtf, features, chroms=None, adjstart=-1, seqprefix=None,
#        fixedwidth=True
#    ):
#    ''' Function to build interval tree covering a specified gtf feature
#    
#    Args:
#        gtf (str)- Path to GTF file.
#        feature (list)- A list of features to find.
#        chroms (dict)- A dictionary of chromosome lengths.
#        adjstart (int)- Adjustment to make to interval start position.
#            A value of -1 will convert 1-based closed intervals to 0-based
#            semi-open intervals.
#        seqprefix (str)- Prefix to add to sequence names.
#    
#    '''
#    # Create output variable
#    intervalDict = {}
#    for feature in features:
#        intervalDict[feature] = find_intervals.FindInterval(chroms=chroms)
#    for entry in gtf_generator(
#            gtf, adjstart=adjstart, seqprefix=seqprefix, fixedwidth=fixedwidth
#        ):
#        if entry.feature != features:
#            continue
#        intervalDict[feature].add_interval(
#            chrom=entry.seqname, start=entry.start, end=entry.end,
#            strand=entry.strand, data=entry.attributes)
#    # Create and return output list
#    outlist = [None] * len(features)
#    for count, feature in enumerate(features):
#        outlist[count] = intervalDict[feature]
#    return(outlist)
