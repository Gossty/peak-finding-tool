import argparse
import pandas as pd
import glob 
from math import log2
from scipy.stats import poisson

WINDOW_LENGTH = 75
GENOME_LENGTH = 2*10**9
THRESHOLD = 1 * 10**(-40)
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
    input_length = gather_data(control_data, args.control)


    get_counts(sample_data, sample_counts)
    get_counts(control_data, control_counts)

    tf_bound = log_fc_filt(sample_counts, control_counts)
    poisson_filt(tf_bound, sample_counts, input_length)

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
        tag_filt = tag[['position', 'read_len', 'strand']]
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
def log_fc_filt(sample_counts, control_counts):
    tf_bound = []
    for index in sample_counts.keys():
        # check if control has this window
        if control_counts.get(index) == None:
            continue
        log_value = log2(sample_counts[index] / control_counts[index])
        if  log_value >= 4:
            tf_bound.append(index)
            # tf_bound[index] = log_value

    return tf_bound


# goes through output from 
def max_fold_filt(log_dict):
    keys = list(log_dict.keys())[0]
    # while
    pass

def poisson_filt(tf_bound, sample_counts, input_length):
    peaks = []
    average_sample_cnt = 0
    average_p_val = 0
    print("total number of tags", input_length)
    print("length of sample counts", len(sample_counts))
    print("length of tf_bound: ", len(tf_bound))
    #lambda for poisson/expected value
    exp = (WINDOW_LENGTH * input_length) / GENOME_LENGTH
    # expected value – shows mean tags for input
    for index in tf_bound:
        average_sample_cnt += sample_counts[index]
        if sample_counts.get(index) == None:
            continue
        p_value = 1 - (poisson.cdf(sample_counts[index], exp))
        average_p_val += p_value
        if p_value < THRESHOLD:
            peaks.append(index)

    average_sample_cnt /= len(tf_bound)
    average_p_val /= len(tf_bound)
    print("average_sample_cnt: ", average_sample_cnt)
    print("average p-value: ", average_p_val)
    print("total number of peaks: ", len(peaks))


    return peaks


main()