from ngs_python.variant import varscan

def geneAnno(
        inFile, outPrefix, annovar, buildver, database
    ):
    '''
    Function perform gene annovar gene annotation on a valid
    annovar input file. Function takes X arguments
    
    '''
    # Build command
    command = [annovar, '-geneanno', '-buildver', buildver, '-outfile',
        outPrefix, inFile, database]
    # Join and return command
    command = ' '.join(command)
    return(command)

def nonSynonomous(
        inFile, outFile, annovar, buildver, database, delete = True
    ):
    # Create file names
    varFile = outFile + '.variant_function'
    exonVarFile = outFile + '.exonic_variant_function'
    logFile = outFile + '.log'
    invalid = outFile + '.invalid_input'
    # Create annovar command
    annovarCommand = geneAnno(inFile = inFile, outPrefix = outFile,
        annovar = annovar, buildver = buildver, database = database)
    # Filter annovar command
    filterCommand = 'awk \'BEGIN{FS="\\t"; OFS="\\t"};{if ($2 ~ '\
        '/^nonsynonymous/) print $9, $3}\' %s > %s' %(
            exonVarFile, outFile)
    # Join commands
    combinedCommand = '%s && %s' %(annovarCommand, filterCommand)
    # Add deletion command if required and return
    if delete:
        combinedCommand += ' && rm %s %s %s %s' %(varFile, exonVarFile,
            logFile, invalid)
    return(combinedCommand)

print nonSynonomous(
    inFile = '/farm/scratch/rs-bio-lif/rabino01/Elza/310123-T1R1.filter.somatic.av',
    outFile = '/farm/scratch/rs-bio-lif/rabino01/Elza/310123-T1R1.filter.somatic.nonsyn',
    annovar = '/farm/babs/redhat6/software/annovar_2015Jun17/annotate_variation.pl',
    buildver = 'mm10',
    database = '/farm/babs/redhat6/software/annovar_2015Jun17/mousedb/',
    delete = True
)
