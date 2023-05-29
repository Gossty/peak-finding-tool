import argparse
import pandas as pd
import glob 
from scipy.stats import poisson
from fuc import pybed
import matplotlib.pyplot as plt
from filters import *
from formating import *

WINDOW_LENGTH = 75
GENOME_LENGTH = 2*10**9
LOCAL_WINDOW = 10000
THRESHOLD = 1 * 10**(-4)
FOLD_VALUE = 4
COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]
COLUMNS_FILT = ['chromosome', 'position', 'read_len', 'strand']

# dictionary stores counts for each window for given starting position
sample_counts = dict()
control_counts = dict()




def main():

    ARGS = arg_parser()
    # formatting 
    formating = Formating(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD)

    # All the dataframes filtered by columns for each chromosome
    sample_df = formating.gather_data( ARGS.tag_directory)
    control_df = formating.gather_data( ARGS.control)

    # Getting all the counts for windows for sample and control
    formating.get_counts(sample_df, sample_counts)
    formating.get_counts(control_df, control_counts)


    # filtering 
    filters = Filters(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                      len(sample_df), len(control_df), FOLD_VALUE)


    # filtering by fold change sample vs control
    tf_bound = filters.fc_filt(sample_counts, control_counts)

    # filtering double counted peaks
    max_filt = filters.max_count_filt(tf_bound, sample_counts)

    # filtering based on the expected number of peaks in control
    poisson_filter = filters.poisson_filt(max_filt, sample_counts)

    # getting the positions of all tags
    dict_tags = formating.get_dict_tags(sample_df)

    # filtering by fold change sample vs LOCAL_WINDOW
    local_filtered = filters.local_filt(sample_counts, poisson_filter, dict_tags)

    # poisson by expected numvber of peaks in LOCAL_WINDOW
    another_poisson = filters.poisson_filt(local_filtered, sample_counts)

    # 
    # # # create a plot to see the p-values
    # # another_poisson.sort()
    # # array_p_value = []
    # # for index in another_poisson: 
    # #     array_p_value.append(sample_counts[index])
    # # exp = (WINDOW_LENGTH * control_length) / GENOME_LENGTH

    # # y = poisson.cdf(array_p_value, mu=exp)
    # # print(max(y))
    # # plt.scatter(array_p_value, y)

        
    # # plt.show()


    # converting to bed file for viewing in IGV
    get_bed(another_poisson)

# getting user input
def arg_parser():
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

    return parser.parse_args()


# building a bed file using pybed library
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