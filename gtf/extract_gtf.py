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
