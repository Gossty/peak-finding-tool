import argparse
import pandas as pd
import glob 
from math import log2
from scipy.stats import poisson
from fuc import pybed

WINDOW_LENGTH = 75
GENOME_LENGTH = 2*10**9
LOCAL_WINDOW = 10000
THRESHOLD = 1 * 10**(-4)
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

    args = parser.parse_args()

    sample_length = gather_data(sample_data, args.tag_directory)
    control_length = gather_data(control_data, args.control)

    df = sample_data['../tests/tagdir/17.tags.tsv']
    dict_tags = get_dict_tags(df)
#
    print("length of sample_length", sample_length)

    get_counts(sample_data, sample_counts)
    get_counts(control_data, control_counts)

    tf_bound = fc_filt(sample_counts, control_counts, sample_length, control_length)
    max_filt = max_count_filt(tf_bound, sample_counts)
    poisson_filter = poisson_filt(max_filt, sample_counts, control_length)

    # local_filtered = local_filt(sample_counts, poisson_filter, dict_tags)

    get_bed(poisson_filter)

# getting all the tag files for tag directory, filtering for necessary columns
# returns total number of tags
def gather_data(data_dict, directory):
    tag_list = []
    tag_list = glob.glob(f"{directory}/*.tsv")
    cnt = 0
    for tag_file in tag_list:
        tag = pd.read_csv(tag_file, sep='\t', header=None)
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


def get_dict_tags(df):    
    pos_list = list(df['position'])

    pos_dict = dict()

    for index in pos_list:
        pos_dict[index] = True
    return pos_dict        


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
def fc_filt(sample_counts, control_counts, sample_length, control_length):
    tf_bound = []
    for index in sample_counts.keys():
        # check if control has this window
        if control_counts.get(index) == None:
            continue
        fold_value = (sample_counts[index] / sample_length) / (control_counts[index] / control_length)
        if  fold_value >= 4:
            tf_bound.append(index)
            # tf_bound[index] = log_value
    tf_bound.sort()
    return tf_bound

# dictionary for helper method
def local_filt(sample_counts, tf_bound, dictionary):
    check = 0
    start_check = 0
    loc_density = 0
    tf_bound_local_filt = []

    for index in tf_bound:

        if check <= 0:
            loc_density = local_density(index, dictionary)
            check = 10000
            start_check = index

        curr_density = sample_counts[index] / WINDOW_LENGTH
        fold_value = curr_density / loc_density
        if  fold_value >= 4:
            tf_bound_local_filt.append(index)
        
        check -= (index - start_check)
        start_check = index
    print("length of local_filt: ", len(tf_bound_local_filt))
    return tf_bound_local_filt


# calculates a tag density for a 10,000 bp window
def local_density(index, dictionary):
    density = 1

    for i in range(index, index + LOCAL_WINDOW + 1):
        if dictionary.get(i) == None:
            continue
        density += 1
    density /= LOCAL_WINDOW
    print("this is local density: ", density)
    return density

# goes through tf_bound and filters to position with max tag count in each window
# window is defined by WINDOW_LENGTH * 2
# returns an array of filtered starting positions
def max_count_filt(tf_bound, sample_counts):
    max_filt = []
    check = tf_bound[0]
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
        # figure out the biggest peak
        max_filt.append(position)
        check += (WINDOW_LENGTH * 2)

    return max_filt

def poisson_filt(tf_bound, sample_counts, control_length):
    peaks = []
    #lambda for poisson/expected value
    exp = (WINDOW_LENGTH * control_length) / GENOME_LENGTH
    # expected value – shows mean tags for input
    for index in tf_bound:
        if sample_counts.get(index) == None:
            continue
        p_value =  1 - poisson.cdf(sample_counts[index], exp)
        if p_value < THRESHOLD:
            peaks.append(index)

    print("total number of peaks in poisson for control: ", len(peaks))



    return peaks

def get_bed(array_filt):

    #removing all values that are not in array_filt
    df = pd.DataFrame()
    chromosome = ["chr17" for i in range(len(array_filt))]
    df.insert(0, "Chromosome", chromosome)
    df.insert(1, "Start", array_filt)

    new_column = [str(int(i) + WINDOW_LENGTH) for i in array_filt]

    #create new column with end position
    df.insert(2, "End", new_column)
    print("length of df", len(df))
    # create a bed file
    bf = pybed.BedFrame.from_frame(meta=[], data=df)
    bf.to_file('example.bed')


main()