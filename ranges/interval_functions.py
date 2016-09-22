import pandas as pd

def merge_overlaps(
        intervals, overlap = 1, return_sorted = True
    ):
    intervals = intervals[['start', 'end']]
    intervals = intervals.sort_values(['end', 'start'])
    intervals['overlap'] = intervals['end'].shift(1) - intervals['start']
    intervals['seperate'] = intervals['overlap'] < overlap
    intervals['group'] = intervals['seperate'].cumsum()
    intervals = intervals.groupby('group')
    # Create and store output dataframe
    output = pd.DataFrame()
    output['start'] = intervals['start'].min()
    output['end'] = intervals['end'].max()
    if return_sorted:
        output = output.sort_values(['start', 'end'])
    output.index = range(output.shape[0])
    return(output)

def merge_overlaps_strand(
        intervals, overlap=1, return_sorted=True, ignore_strand=False
    ):
    # Merge intervals while ignoring strand
    if ignore_strand:
        intervals = intervals[['start', 'end']]
        intervals = merge_overlaps(
            intervals, overlap, False
        )
        intervals['strand'] = '*'
    # Merge intervals while including strand
    else:
        intervals = intervals[['start', 'end', 'strand']]
        intervals = intervals.groupby('strand')
        intervals = intervals.apply(
            lambda x: merge_overlaps(x, overlap, False)
        )
        intervals = intervals.reset_index('strand')
        intervals = intervals[['start', 'end', 'strand']]
    # Process and return output intervals
    if return_sorted:
        intervals = intervals.sort_values(['start', 'end'])
    intervals.index = range(intervals.shape[0])
    return(intervals)

def merge_overlaps_chrom(
    intervals, overlap=1, return_sorted=True, ignore_strand=False
):
    # Generate intervals data frame and split on strand
    if ignore_strand:
        intervals = intervals[['chr', 'start', 'end']]
        intervals['strand'] = '*'
    else:
        intervals = intervals[['chr', 'start', 'end', 'strand']]
    # Group intervals by strand and merge overlaps
    intervals = intervals.groupby('chr')
    intervals = intervals.apply(
        lambda x: merge_overlaps_strand(x, overlap, False)
    )
    intervals = intervals.reset_index('chr')
    # Process and return intervals
    print(intervals)
    if return_sorted:
        intervals = intervals.sort_values(['chr', 'start', 'end', 'strand'])
    intervals.index = range(intervals.shape[0])
    return(intervals)

x = pd.DataFrame()
x['chr'] = ['chr1', 'chr1', 'chr2', 'chr2', 'chr3']
x['start'] = [0, 5, 0, 5, 0]
x['end'] = [10, 15, 10, 15, 10]
x['strand'] = ['+', '-', '-', '-', '*']
print(x)
y = merge_overlaps_chrom(x, ignore_strand=True)
print(y)
