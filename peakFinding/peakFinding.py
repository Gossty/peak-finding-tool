import argparse
import pandas as pd
import glob 
from math import log2
from scipy.stats import poisson
from fuc import pybed

WINDOW_LENGTH = 75
GENOME_LENGTH = 2*10**9
THRESHOLD = 1 * 10**(-60)
COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]

# dictionary stores counts for each window for given starting position
sample_counts = dict()
control_counts = dict()

# All the filtered dataframes for each chromosome
sample_data = dict()
control_data = dict()



def main():
    parser = argparse.ArgumentParser(
        prog='peakFinding',
        description='Finds Peaks in the provided tag directory'
    )
    parser.add_argument('tag_directory', help="tag directory for analysis")
    parser.add_argument("control", help="control of peak finding")
    parser.add_argument('-style', help='<factor> for TFs; <histone> for histone modifications')
    parser.add_argument('-o',help='output directory')
    parser.add_argument('-poisson',help='manually set the threshold for poisson')
    parser.add_argument('-fold', help='manually set the fold change for peak detection')
    parser.add_argument('-fragLen',help='manually specify the length of window')
    parser.add_argument('-L', help='manually set the scope for local filtering')

    args = parser.parse_args()

    sample_length = gather_data(sample_data, args.tag_directory)
    input_length = gather_data(control_data, args.control)


    get_counts(sample_data, sample_counts)
    get_counts(control_data, control_counts)

    tf_bound = fc_filt(sample_counts, control_counts)
    tf_bound_filt = max_count_filt(tf_bound, sample_counts)
    poisson_filt(tf_bound_filt, sample_counts, input_length)

    get_bed(sample_data, tf_bound_filt)

# getting all the tag files for tag directory, filtering for necessary columns
# returns total number of tags
def gather_data(data_dict, directory):
    tag_list = []
    tag_list = glob.glob(f"{directory}/*.tsv")
    cnt = 0
    for tag_file in tag_list:
        tag = pd.read_csv(tag_file, sep='\t')
        tag.columns = COLUMNS
        # removing unnecesary
        tag_filt = tag[['chromosome', 'position', 'read_len', 'strand']]
        data_dict[tag_file] = tag_filt
        cnt += len(tag_filt)
    return cnt


# runs through all the tags for each sample
# for each read in the sample updates counts based on overlap
def get_counts(data, dictionary):
    for tag_file, tag_filt in data.items():
        for index, row in tag_filt.iterrows():
            overlap(dictionary, row['position'], row['read_len'], row['strand'])


# overlap – increments all the values in the dictionary of windows 
# tag – start index of tag, tag_len – length of tag, dictionary – read
def overlap(dictionary, tag, tag_len, strand):
        center = int(tag + tag_len / 2)
        # the range of for loop is representing the starting positions of the window
        if strand: # checks for strand direction
            for x in range(center - WINDOW_LENGTH, center + 1):
                # checking if the position exists; if yes – increment
                if dictionary.get(x) == None:
                    dictionary[x] = 1
                else:
                    dictionary[x] += 1
        else:
            # the range of for loop is representing the starting positions of the window
            # range is reversed
            for x in range(center + WINDOW_LENGTH, center + 1, -1):
                # checking if the position exists; if yes – increment
                if dictionary.get(x) == None:
                    dictionary[x] = 1
                else:
                    dictionary[x] += 1


# filtering based of log fold change.
# tf_bound keeps track of positions where TFs bound based on whether log2 fold change < 4
def fc_filt(sample_counts, control_counts):
    tf_bound = []
    for index in sample_counts.keys():
        # check if control has this window
        if control_counts.get(index) == None:
            continue
        fold_value = sample_counts[index] / control_counts[index]
        if  fold_value >= 4:
            tf_bound.append(index)
            # tf_bound[index] = log_value
    tf_bound.sort()
    return tf_bound


# goes through tf_bound and filters to position with max tag count in each window
# window is defined by WINDOW_LENGTH * 2
# returns an array of filtered starting positions
def max_count_filt(tf_bound, sample_counts):
    tf_bound_filt = []
    check = tf_bound[0]
    yang_peak = tf_bound[0]
    biggest_peak = -1
    # iterating through starting position of tf_bound
    for index in tf_bound:
        # continue iterating until the check
        if index < check:
            continue
        check = index

        start = index
        position = start
        max_dict = dict()
        for i in range(start, start + (WINDOW_LENGTH * 2) + 1):
            if sample_counts.get(i) == None:
                continue   
            max_dict[i] = sample_counts.get(i)
        
        position = max(max_dict)
        local_max = max(max_dict.values())
        # figure out the biggest peak
        if biggest_peak < local_max:
            biggest_peak = local_max
            yang_peak = position
        tf_bound_filt.append(position)
        check += (WINDOW_LENGTH * 2)
    print("filtering tf_bound by maximum value ",len(tf_bound_filt))
    print("THE YANG PEAK POGGERS 🥰", yang_peak)
    print("MAXIMUM TAG COUNT 🔝", biggest_peak)
    return tf_bound_filt

def poisson_filt(tf_bound, sample_counts, input_length):
    peaks = []
    average_sample_cnt = 0
    average_p_val = 0
    #lambda for poisson/expected value
    exp = (WINDOW_LENGTH * input_length) / GENOME_LENGTH
    # expected value – shows mean tags for input
    for index in tf_bound:
        average_sample_cnt += sample_counts[index]
        if sample_counts.get(index) == None:
            continue
        p_value =  1 - poisson.cdf(sample_counts[index], exp)
        average_p_val += p_value
        if p_value < THRESHOLD:
            peaks.append(index)

    average_sample_cnt /= len(tf_bound)
    average_p_val /= len(tf_bound)
    print("average_sample_cnt: ", average_sample_cnt)
    print("average p-value: ", average_p_val)
    print("total number of peaks: ", len(peaks))
    print(peaks[:100])



    return peaks

def get_bed(df, array_filt):
    # df[df['A'].isin([3, 6])]
    
    df_filt = df['../tests/tagdir/17.tags.tsv']
    #removing all values that are not in array_filt

    print("before is in: ", len(df_filt))

    df_filt = df_filt[df_filt['position'].isin(array_filt)]
    df_filt = df_filt[['chromosome', 'position']]

    #changing the chromosome name
    df_filt.loc[df_filt['chromosome'] == 17, 'chromosome'] = "chr17"

    startings = list(df_filt['position'])

    #create new column with end position
    new_column = [str(int(i) + WINDOW_LENGTH) for i in startings]
    df_filt.insert(2, "end", new_column)
    print("after is in method: ", len(df_filt))
    df_filt.columns = ['Chromosome', 'Start', 'End']
    # create a bed file
    bf = pybed.BedFrame.from_frame(meta=[], data=df_filt)
    bf.to_file('example.bed')


main()