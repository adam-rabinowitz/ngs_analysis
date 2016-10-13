import re
import pandas as pd

def extract_exons_bed(gene, gtf, bed, identifier='id', source=None):
    # Compile regular expression for identifer
    if identifier == 'id':
        pattern = 'gene_id "{}"'.format(gene)
    elif identifier == 'name':
        pattern = 'gene_name "{}"'.format(gene)
    else:
        raise ValueError('Identifier must be id or name')
    idRegx = re.compile(pattern)
    # Compile regular expression for exon number
    noRegx = re.compile('exon_number "(\\d+)"')
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
            exonSearch = noRegx.search(lineData[8])
            if exonSearch:
                name = '{}_exon_{}_{}'.format(
                    gene,
                    exonSearch.group(1),
                    lineData[1]
                )
            else:
                name = '{}_exon_{}'.format(
                    gene,
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
    ''' Function to generate bed file listing the first transcribed
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
            # Extract gene id
            geneID = idRegx.search(lineData[8]).group(1)
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
                geneID,
                lineData[6]
            )
            outFile.write(outLine)

extract_tss(
    '/farm/scratch/rs-bio-lif/rabino01/Ascl1/genomeData/Mus_musculus.GRCm38.80.exon.gtf',
    '/farm/scratch/rs-bio-lif/rabino01/Ascl1/genomeData/Mus_musculus.GRCm38.80.tss.bed'
)
