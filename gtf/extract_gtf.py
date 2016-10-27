import re
import pandas as pd

def extract_exons_bed(
        ident, ident_type, gtf, bed, source=None
    ):
    # Compile regular expression for identifer
    if ident_type == 'gene':
        pattern = 'gene_id "{}"'.format(ident)
    elif ident_type == 'tran':
        pattern = 'transcript_id "{}"'.format(ident)
    elif ident_type == 'name':
        pattern = 'gene_name "{}"'.format(ident)
    else:
        raise ValueError('ident_type must be one of: gene, tran or name')
    idRegx = re.compile(pattern)
    # Compile regular expression for exon number
    exonRegx = re.compile('exon_number "(\\d+)"')
    tranRegx = re.compile('transcript_id "(ENSMUST\\d+)"')
    # Loop through output
    count = 0
    with open(gtf) as inFile, open(bed, 'w') as outFile:
        for line in inFile:
            if line.startswith('#'):
                continue
            lineData = line.strip().split('\t')
            if lineData[2] != 'exon':
                continue
            if source and lineData[1] != source:
                continue
            if not idRegx.search(lineData[8]):
                continue
            name = '{}_exon_{}_{}'.format(
                tranRegx.search(lineData[8]).group(1),
                exonRegx.search(lineData[8]).group(1),
                lineData[1]
            )
            bedLine = '{}\t{}\t{}\t{}\t0\t{}\n'.format(
                lineData[0],
                int(lineData[3]) - 1,
                lineData[4],
                name,
                lineData[6]
            )
            print(bedLine)
            outFile.write(bedLine)

def extract_tss(gtf, bed, source=None):
    ''' Function to identrate bed file listing the first transcribed
    base of all transcripts in a bed file.
    
    Args:
        gtf (str)- Full path to input gtf file.
        bed (str)- Full path to output bed file.
        source (str)- Value to match in source column of gtf.
    
    '''
    noRegx = re.compile('exon_number "(\\d+)"')
    idRegx = re.compile('transcript_id "(\\w+)"')
    with open(gtf) as inFile, open(bed, 'w') as outFile:
        for line in inFile:
            # Skip header lines
            if line.startswith('#'):
                continue
            # Extract line data and check source
            lineData = line.strip().split('\t')
            if source and lineData[1] != source:
                continue
            # Extract and check exon number
            exonNumber = noRegx.search(lineData[8]).group(1)
            if exonNumber != '1':
                continue
            # Extract ident id
            identID = idRegx.search(lineData[8]).group(1)
            # Find exon start
            if lineData[6] == '+':
                tss = int(lineData[3])
            elif lineData[6] == '-':
                tss = int(lineData[4])
            else:
                raise ValueError("Strand must be '+' or '-'")
            # Create output line
            start = tss - 1
            end = tss
            outLine = '{}\t{}\t{}\t{}\t0\t{}\n'.format(
                lineData[0],
                start,
                end,
                identID,
                lineData[6]
            )
            outFile.write(outLine)
