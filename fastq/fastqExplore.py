import collections
import gzip
import sys
import argparse
from fuzzywuzzy import fuzz

def parseFasta(fasta, length):
    ''' Function generates a dictionary from a FASTA file where the keys
    are the sequence names and the values are the sequences.
    
    Args:
        fasta (str)- Path to the fasta file.
        length (int)- Length of the sequences to extract.
    
    Returns:
        fastaDict (dict)- Keys are the sequence names and values are the
            sequences.
    
    '''
    # Create valriables
    fastaDict = {}
    name = None
    # Open file and loop through entries
    with open(fasta) as infile:
        for line in infile:
            # Extract sequence names and store previous entries
            if line.startswith('>'):
                if name is not None:
                    fastaDict[name] = sequence
                name = line[1:].strip()
                sequence = ''
            # Append sequence
            else:
                line = line.strip('^')
                sequence = sequence + line.strip()
        # Add final sequence and return dictionary
        if name is not None:
            fastaDict[name] = sequence
    return(fastaDict)

def findBestMatch(query, sequenceDict, trim=True):
    ''' Find best match between a query sequence and sequence in a
    dictionary.
    
    Args:
        query (str)- Query sequence.
        sequenceDict (dict)- Dictionary where values are the the sequences
            against which to match the sequence.
        trim (bool)- If True then the longest sequence is trimmed to
            the length of the shortest sequence prior to matching.
    
    Returns:
        bestMatch (list)- List of the keys of the best matching sequences.
        bestScore (float)- Score of the best match.
    
    '''
    # Create output variables
    bestMatch = []
    bestScore = 0
    # Loop through
    for name, sequence in sequenceDict.items():
        # Trim sequences, if required, and generate score
        if trim:
            length = min(len(query), len(sequence))
            score = fuzz.ratio(query[:length], sequence[:length])
        else:
            score = fuzz.ratio(query, sequence)
        # Store best score
        if score > bestScore:
            bestScore = score
            bestMatch = [name]
        # Else append equal best score
        elif score == bestScore:
            bestMatch.append(name)
    # Return names of nest matches and score
    return((bestMatch, bestScore))

def fastqSeqIterator(
        fastq, number, length
    ):
    ''' Generator returning the desired number of sequences froma FASTQ
    file. Function presumes that each entry in the FASTQ file is four
    lines.
    
    Args:
        fastq - Path to FASTQ file.
        number - Number of sequences to return.
        length - Length of the sequences to return.
    
    Yields:
        startseq - Start sequence of each entry.
    
    '''
    # Create open function
    if fastq.endswith('.gz'):
        openfunc = gzip.open
    else:
        openfunc = open
    # Open fastq file and loop through lines
    with openfunc(fastq) as infile:
        for count, line in enumerate(infile):
            if (count % 4) == 1:
                startseq = line[:length]
                yield(startseq)
            if (count / 4) >= number:
                break

def findCommonFastqStartSequences(
        fastq, number=100000, length=10, report=10, fasta=None, trim=True
    ):
    ''' Function to identify, and annotate, common sequences at the start
    of reads within a FASTQ file.
    
    Args:
        fastq (str)- Path to fastq file.
        number (int)- Number of FASTQ reads to parse.
        length (int)- Length of sequence to report.
        report (int)- Number of sequences to report.
        fasta (str)- Path to FASTA file contianing sequences against which
            to match the sequences.
        trim (bool)- Whether to trim query and refence sequences to shortest
            common length prior to matching.
    
    '''
    # Parse fasta dictionary
    if fasta is not None:
        fastaDict = parseFasta(
            fasta=fasta, length=length)
    else:
        fastaDict = None
    # Extract start sequences
    startSeq = collections.Counter()
    for seq in fastqSeqIterator(
            fastq=fastq, number=number, length=length
        ):
        startSeq[seq] += 1
    # Create outut string
    if fastaDict is None:
        outstring = '{}\t{}\t{:.4f}'
    else:
        outstring = '{}\t{}\t{:.4f}\t{}\t{}'
    # Create variables for tallying sequences
    total = sum(startSeq.values())
    reportedCount = 0
    reportedRatio = 0.0
    # Loop through common sequences and report
    for sequence, count in startSeq.most_common(report):
        # Calculate ratio and store progressive counts
        ratio = float(count) / total
        reportedCount += count
        reportedRatio += ratio
        # Report results, with or without fasta matching
        if fastaDict is None:
            print(outstring.format(sequence, count, ratio))
        else:
            match, score = findBestMatch(
                query=sequence, sequenceDict=fastaDict, trim=trim)
            match = ','.join(match)
            print(outstring.format(sequence, count, ratio, match, score))
    # Report remianing sequences
    otherCount = total - reportedCount
    otherRatio = 1 - reportedRatio
    if fastaDict is None:
        print(outstring.format('other', otherCount, otherRatio))
    else:
        print(outstring.format('other', otherCount, otherRatio, '', ''))
