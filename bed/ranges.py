def independent(df):
    # Extract data for each chromosome
    for c in np.unique('chr'):
        # Build an array on interleaved start and stop values
        chrData = df[df['chr'] == c]
        ranges = np.ravel(chrData[['start','end']])
        # Return false if ranges are not ordered and non-overlapping
        diff = ranges[1:] - ranges[:-1]
        if not all(diff > 0):
            return(False)
    # Return true if all chromosomes are ordered
    return(True)

def overlap(query, target):
    for 
