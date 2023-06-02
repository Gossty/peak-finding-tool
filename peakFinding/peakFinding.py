import argparse
import pandas as pd
import matplotlib.pyplot as plt
from filters import *
from formating import *
import sys


# dictionary stores counts for each window for given starting position
sample_counts = dict()
control_counts = dict()

# getting user input
def arg_parser():
    parser = argparse.ArgumentParser(
        prog='peakFinding',
        description='Finds Peaks in the provided tag directory'
    )

    parser.add_argument('tag_directory', help="tag directory for analysis")
    parser.add_argument("-control", help="control of peak finding")
    parser.add_argument('-style', help='<factor> for TFs; <histone> for histone modifications')
    parser.add_argument('-o',help='output directory')
    parser.add_argument('-poisson',help='manually set the threshold for poisson')
    parser.add_argument('-fold', help='manually set the fold change for peak detection')
    parser.add_argument('-fragLen',help='manually specify the length of window')
    parser.add_argument('-L', help='manually set the scope for local filtering')

    return parser.parse_args()


ARGS = arg_parser()

WINDOW_LENGTH = 75 if ARGS.fragLen == None else ARGS.fragLen
GENOME_LENGTH = 2*10**9
LOCAL_WINDOW = 10000 if ARGS.L == None else ARGS.L
THRESHOLD = 1 * 10**(-4) if ARGS.poisson == None else ARGS.poisson
FOLD_VALUE = 4 if ARGS.fold == None else ARGS.fold
COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]
COLUMNS_FILT = ['chromosome', 'position', 'read_len', 'strand']
OUT_FILE='peak.txt' if ARGS.o == None else ARGS.o




def main():
    flag= 0
    #Updating the values based user input
    if hasattr(ARGS,"tag_directory")==False:
        raise Exception("Please provide tag_directory. Type 'peakFinding -help' for the usage of this program")
    
    flag = 0 if ARGS.control == None else 1
    
    #-style is not useful at this stage since histone mode is not implemented
    
    # formatting 
    formating = Formating(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD)

    # All the dataframes filtered by columns for each chromosome
    sample_df = formating.gather_data( ARGS.tag_directory)

    # Getting all the counts for windows for sample and control
    formating.get_counts(sample_df, sample_counts)

    sorted_positions = list(sample_counts.keys())
    sorted_positions.sort()

    if flag==0:
        # filtering 
        filters = Filters(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                      len(sample_df), len(sample_df), FOLD_VALUE)

        # filtering double counted peaks
        peaks_output = filters.max_count_filt(sorted_positions, sample_counts)
        array_for_graph = peaks_output
        PUTATIVE_PEAKS = len(peaks_output)

        # getting the positions of all tags
        dict_tags = formating.get_dict_tags(sample_df)

        # filtering by fold change sample vs LOCAL_WINDOW
        peaks_output = filters.local_filt(sample_counts, peaks_output, dict_tags)
        PUTATIVE_BY_LOC = len(peaks_output)

        peaks_output = filters.poisson_filt(peaks_output, sample_counts)
        peaks_output.sort()
        peak_stats(peaks_output, sample_counts, sample_df, ARGS, 
                   PUTATIVE_PEAKS, PUTATIVE_BY_LOC, OUT_FILE)


    else:
        
        # All the dataframes filtered by columns for each chromosome
        control_df = formating.gather_data( ARGS.control)

        # Getting all the counts for windows for sample and control
        formating.get_counts(control_df, control_counts)

        filters = Filters(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                      len(sample_df), len(control_df), FOLD_VALUE)

        # filtering double counted peaks
        peaks_output = filters.max_count_filt(sorted_positions, sample_counts)
        array_for_graph = peaks_output
        PUTATIVE_PEAKS = len(peaks_output)

        # filtering by fold change sample vs control
        peaks_output = filters.fc_filt(sample_counts, control_counts, peaks_output)
        PUTATIVE_BY_INPUT = len(peaks_output)

        # poisson by expected numvber of peaks in LOCAL_WINDOW
        peaks_output = filters.poisson_filt(peaks_output, sample_counts)

        # getting the positions of all tags
        dict_tags = formating.get_dict_tags(sample_df)

        # filtering by fold change sample vs LOCAL_WINDOW
        peaks_output = filters.local_filt(sample_counts, peaks_output, dict_tags)
        PUTATIVE_BY_LOC = len(peaks_output)

        # filtering based on the expected number of peaks in control
        peaks_output = filters.poisson_filt(peaks_output, sample_counts)
        peaks_output.sort()
        peak_stats(peaks_output, sample_counts, sample_df, ARGS, 
                   PUTATIVE_PEAKS, PUTATIVE_BY_LOC, OUT_FILE, PUTATIVE_BY_INPUT,sample_df)



    # false_peaks(sample_counts, control_counts, 34000123, WINDOW_LENGTH)

    

    # x = list(range(0, 151))
    # bar_chart = [0] * 151
    # for index in array_for_graph:
    #     bar_chart[sample_counts[index]] +=1
    # plt.bar(x, bar_chart)
    # plt.show()


    # converting to bed file for viewing in IGV
    get_bed(peaks_output)


def false_peaks(sample_counts, control_counts, position, WINDOW_LENGTH):
    cnt = 0
    for index in range(position, position + WINDOW_LENGTH * 2 + 1):
        if sample_counts.get(index) == None:
            cnt += 1
        print(f"""position: {index}, 
        sample: {sample_counts[index]}, control: {control_counts[index]}, 
        fold_change: {sample_counts[index] / control_counts[index]}""")

    if cnt >= WINDOW_LENGTH * 2:
        print("No windows?")

def peak_stats(peaks_output, sample_counts, sample_df, input,
               putative_peaks, putative_by_loc, file_path, putative_by_input=-1, control_df=-1):

    
    total_peaks = len(peaks_output)
    peak_size = WINDOW_LENGTH
    minimum_distance_peaks = GENOME_LENGTH
    genome_size = GENOME_LENGTH
    fold_val = FOLD_VALUE
    p_value_req = THRESHOLD
    local_wind = LOCAL_WINDOW

    total_tags = len(sample_df)
    tags_in_peaks = 0
    tags_per_bp = 0

    if type(control_df) == 'DataFrame':
        tags_per_bp = len(control_df) / GENOME_LENGTH
    else:
        tags_per_bp = len(sample_df) / GENOME_LENGTH
    
    exp_tags_per_peak = WINDOW_LENGTH * tags_per_bp


    tags_for_normalization = len(sample_df)
    command = ' '.join(sys.argv)


    prev_peak = peaks_output[0]
    for index in peaks_output:
        tags_in_peaks += sample_counts[index]
        delta = index - prev_peak

        if delta != 0 and delta < minimum_distance_peaks:
            minimum_distance_peaks = delta

        prev_peak = index

    ip_efficiency = (tags_in_peaks / total_tags) * 100

    try:
        file = open(file_path, 'w')
    except FileNotFoundError:
        print("The given output path is invalid")

    # Write content to the file
    file.write('# YSY Peaks \n')
    file.write('# Peak finding parameters: \n')
    file.write(f'# tag directory = {input.tag_directory} \n')
    file.write(f'# total peaks = {total_peaks}\n')
    file.write(f'# peak size = {peak_size}\n')
    file.write(f'# peaks found using tags on both strands \n')
    file.write(f'# minimum distance between peaks = {minimum_distance_peaks}\n')
    file.write(f'# genome size = {genome_size} \n')
    file.write(f'# total tags = {total_tags} \n')
    file.write(f'# total tags in peaks = {tags_in_peaks}\n')
    file.write(f'# approximate IP efficiency = {ip_efficiency }% \n')
    file.write(f'# tags per bp = {tags_per_bp}\n')
    file.write(f'# expected tags per peak = {exp_tags_per_peak}\n')
    file.write(f'# effective number of tags used for normalization = {tags_for_normalization}\n')
    file.write(f'# number of putative peaks = {putative_peaks} \n')
    file.write(f'#  \n')

    if type(control_df) != 'DataFrame':
        file.write(f'# input tag directory = {input.control} \n')
        file.write(f'# Fold over input required = {fold_val}\n')
        file.write(f'# Poisson p-value over input required = {p_value_req} \n')
        file.write(f'# Putative peaks filtered by input = {putative_by_input}  \n')
        file.write(f'#  \n')
    
    file.write(f'# size of region used for local filtering = {local_wind} \n')
    file.write(f'# Fold over local region required = {fold_val} \n')
    file.write(f'# Poisson p-value over local region required = {p_value_req} \n')
    file.write(f'# Putative peaks filtered by local filter = {putative_by_loc}  \n')
    file.write(f'#  \n')

    file.write(f'# cmd = {command} \n')

    # Close the file
    file.close()




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

    # create a bed file
    file = open(f'output.bed', 'w')
    with file as f:
        df.apply(write_row, args=(file,), axis=1)

def write_row(row, file):
    
    file.write(f"{row['Chromosome']}\t{row['Start']}\t{row['End']}")
    file.write('\n')



main()