import argparse
import pandas as pd
import os
from .filters import *
from .formating import *
import sys

# constants for printing
RED = '\033[91m'
RESET = '\033[0m'

# getting user input
def arg_parser():
    parser = argparse.ArgumentParser(
        prog='peakFinding',
        description='Find peaks in a provided tag directory. It takes as input a tag directory containing tag files, and it optionally takes a control directory for peak finding as well as other parameters. The code performs various filtering and analysis steps to identify peaks and outputs the results in a BED file format and a txt file with statistical information.'
    )

    parser.add_argument('tag_directory', help="tag directory for analysis")
    parser.add_argument("-control", help="control directory for peak finding")
    parser.add_argument('-o',help='output directory')
    parser.add_argument('-poisson',help='manually set the threshold for poisson')
    parser.add_argument('-fold', help='manually set the fold change for peak detection')
    parser.add_argument('-fragLen',help='manually specify the length of window')
    parser.add_argument('-L', help='manually set the scope for local filtering')
    return parser.parse_args()

# assigning user input
ARGS = arg_parser()
WINDOW_LENGTH = 75 if ARGS.fragLen == None else int(ARGS.fragLen)
GENOME_LENGTH = 2*10**9
LOCAL_WINDOW = 10000 if ARGS.L == None else int(ARGS.L)
THRESHOLD = 1 * 10**(-4) if ARGS.poisson == None else float(ARGS.poisson) 
FOLD_VALUE = 4 if ARGS.fold == None else float(ARGS.fold)
OUT_DIRECTORY='./' if ARGS.o == None else str(ARGS.o)
COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]
COLUMNS_FILT = ['chromosome', 'position', 'read_len', 'strand']

if os.path.exists(OUT_DIRECTORY) == False:
    raise Exception(RED + "Error: Path to output directory doesn't exist" + RESET)


# main function of the program that runs the filters according to the user input
# outputs peaks.txt and peaks.bed files
def main():
    if ARGS.tag_directory == None:
        raise Exception(RED + "Error: Please provide tag_directory. Type 'peakFinding -help' for the usage of this program" + RESET)
    if os.path.exists(ARGS.tag_directory) == False:
        raise Exception(RED + "Error: Path to tag_directory doesn't exist" + RESET)
    
    FLAG = 0 if ARGS.control == None else 1
    
    
    # formatting 
    formating = Formating(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD)

    # All the dataframes filtered by columns for each chromosome
    total_dict_sample = dict()
    sample_df = formating.gather_data(ARGS.tag_directory, total_dict_sample)


    total_output = dict()

    # counters for statistics in peaks.txt
    putative_peaks = 0
    putative_by_input = 0
    putative_by_loc = 0
    tags_in_peaks = 0
    
    if FLAG == 0:
        for chromosome in total_dict_sample.keys():
            # Getting dataframe with tags for the current chromosome
            chr_df = sample_df.loc[(sample_df['chromosome'] == chromosome)]
            sample_counts = dict()
            # Getting all the counts for windows for sample and control
            formating.get_counts(chr_df, sample_counts)
            print("filtering chromosome:", chromosome)
            sorted_positions = list(sample_counts.keys())
            # filtering 
            filters = Filters(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                        len(sample_df), len(sample_df), FOLD_VALUE)

            # filtering double counted peaks
            peaks_output = filters.max_count_filt(sorted_positions, sample_counts)
            putative_peaks += len(peaks_output)

            # getting the positions of all tags
            dict_tags = formating.get_dict_tags(chr_df)

            # filtering by fold change sample vs LOCAL_WINDOW
            peaks_output = filters.local_filt(sample_counts, peaks_output, dict_tags)
            putative_by_loc += len(peaks_output)

            peaks_output, tags_in_peaks = filters.poisson_filt(peaks_output, sample_counts, tags_in_peaks)

            peaks_output.sort()
            total_output[chromosome] = peaks_output

        # printing txt file with all the stats
        peak_stats(total_output, tags_in_peaks, sample_df, ARGS, 
                putative_peaks, putative_by_loc, OUT_DIRECTORY)

    else:

        if os.path.exists(ARGS.control) == False:
            raise Exception(RED + "Error: Path to control_directory doesn't exist" + RESET)
        
        total_dict_input = dict()
        control_df = formating.gather_data(ARGS.control, total_dict_input)
        
        for chromosome in total_dict_sample.keys():
            chr_df = sample_df.loc[(sample_df['chromosome'] == chromosome)]
            sample_counts = dict()
            formating.get_counts(chr_df, sample_counts)
            
            if total_dict_input.get(chromosome) == None:
                continue

            control_chr_df = control_df.loc[(control_df['chromosome'] == chromosome)]
            control_counts= dict()
            # All the dataframes filtered by columns for each chromosome
            formating.get_counts(control_chr_df, control_counts)


            print("filtering chromosome:", chromosome)

            sorted_positions = list(sample_counts.keys())
            # Getting all the counts for windows for sample and control
            filters = Filters(WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                        len(sample_df), len(control_df), FOLD_VALUE)
            

            # filtering double counted peaks
            peaks_output = filters.max_count_filt(sorted_positions, sample_counts)
            putative_peaks += len(peaks_output)

            # filtering by fold change sample vs control
            peaks_output = filters.fc_filt(sample_counts, control_counts, peaks_output)
            putative_by_input += len(peaks_output)

            # poisson by expected numvber of peaks in LOCAL_WINDOW
            peaks_output = filters.poisson_filt(peaks_output, sample_counts)

            # getting the positions of all tags
            dict_tags = formating.get_dict_tags(chr_df)

            # filtering by fold change sample vs LOCAL_WINDOW
            peaks_output = filters.local_filt(sample_counts, peaks_output, dict_tags)
            putative_by_loc += len(peaks_output)

            # filtering based on the expected number of peaks in control
            peaks_output, tags_in_peaks = filters.poisson_filt(peaks_output, 
                                                  sample_counts, tags_in_peaks)
            peaks_output.sort()
            total_output[chromosome] = peaks_output
        
        # printing txt file with all the stats
        peak_stats(total_output, tags_in_peaks, sample_df, ARGS, 
                putative_peaks, putative_by_loc, OUT_DIRECTORY, putative_by_input,sample_df)



    



    # converting to bed file for viewing in IGV
    get_bed(total_output, OUT_DIRECTORY)

    sys.exit(0)


# method that returns peaks.txt file through the statistical data
def peak_stats(total_output, tags_in_peaks, sample_df, input, putative_peaks, 
                putative_by_loc, file_path, putative_by_input=-1, control_df=-1):

    
    total_peaks = 0
    peak_size = WINDOW_LENGTH
    minimum_distance_peaks = GENOME_LENGTH
    genome_size = GENOME_LENGTH
    fold_val = FOLD_VALUE
    p_value_req = THRESHOLD
    local_wind = LOCAL_WINDOW

    total_tags = len(sample_df)
    tags_per_bp = 0
    for chromosome in total_output.keys():
        peaks_output = total_output[chromosome]
        total_peaks += len(peaks_output)
        if len(peaks_output) == 0:
            continue
        prev_peak = peaks_output[0]

        # calculate minimum distance between peaks
        for index in peaks_output:
            delta = index - prev_peak

            if delta != 0 and delta < minimum_distance_peaks:
                minimum_distance_peaks = delta
            prev_peak = index

    if type(control_df) == 'DataFrame':
        tags_per_bp = len(control_df) / GENOME_LENGTH
    else:
        tags_per_bp = len(sample_df) / GENOME_LENGTH
    
    exp_tags_per_peak = WINDOW_LENGTH * tags_per_bp


    tags_for_normalization = len(sample_df)
    command = ' '.join(sys.argv)


    ip_efficiency = (tags_in_peaks / total_tags) * 100

    try:
        file = open(f'{file_path}/peaks.txt', 'w')
    except FileNotFoundError:
        print(RED+"The given output path is invalid"+RESET)

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
    print(f"peaks.txt output goes to this directory: {file_path}")
    # Close the file
    file.close()




# building a bed file using pybed library
def get_bed(total_output, OUT_DIRECTORY):
    file = open(f'{OUT_DIRECTORY}/peaks.bed', 'w')

    for chromosome in total_output.keys():
        file = open(f'{OUT_DIRECTORY}/peaks.bed', 'a')

        array_filt = total_output[chromosome]
        if len(array_filt) == 0:
            continue


        #removing all values that are not in array_filt
        df = pd.DataFrame()
        chromosome_arr = []
        if 'chr' in str(chromosome):
            chromosome_arr = [f"{chromosome}" for i in range(len(array_filt))]
        else:
            chromosome_arr = [f"chr{chromosome}" for i in range(len(array_filt))]
        df.insert(0, "Chromosome", chromosome_arr)
        df.insert(1, "Start", array_filt)

        new_column = [str(int(i) + WINDOW_LENGTH) for i in array_filt]

        #create new column with end position
        df.insert(2, "End", new_column)

        # create a bed file
        with file as f:
            df.apply(write_row, args=(file,), axis=1)

    print(f"peaks.bed output goes to this directory: {OUT_DIRECTORY}")

# helper method for get_bed
def write_row(row, file):
    file.write(f"{row['Chromosome']}\t{row['Start']}\t{row['End']}")
    file.write('\n')

if __name__ == "__main__":
    main()
