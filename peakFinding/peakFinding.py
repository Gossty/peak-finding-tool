import argparse
import pandas as pd
import glob 
from math import log2

WINDOW_LENGTH = 75
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

    args = parser.parse_args()

    gather_data(sample_data, args.tag_directory)
    gather_data(control_data, args.control)

    get_counts(sample_data, sample_counts)
    get_counts(control_data, control_counts)

    log_fc_filt(sample_counts, control_counts)


# getting all the tag files for tag directory, filtering for necessary columns
def gather_data(data_dict, directory):
    tag_list = []
    tag_list = glob.glob(f"{directory}/*.tsv")
    for tag_file in tag_list:
        tag = pd.read_csv(tag_file, sep='\t')
        tag.columns = COLUMNS
        # removing unnecesary
        tag_filt = tag[['position', 'read_len', 'strand']]
        # ensure that it is sorted by position
        # tag_filt.sort_values(by="position", ascending=True, inplace=True)
        data_dict[tag_file] = tag_filt



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


def log_fc_filt(sample_counts, control_counts):
    log_dict = dict()
    cnt = 0
    print("length of sample_counts before filtering", len(sample_counts))
    for index in sample_counts.keys():
        if control_counts.get(index) == None:
            sample_counts.pop(index)
            continue
        elif log2(sample_counts[index] / control_counts[index]) < 4:
            cnt += 1
            sample_counts.pop(index)
            control_counts.pop(index)
        
    print("less than 4: ", cnt)

    print("length of sample_counts after filtering", len(sample_counts))


    # return log_dict


main()