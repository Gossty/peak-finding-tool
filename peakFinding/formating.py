import glob 
import pandas as pd
from fuc import pybed



class Formating:
    def __init__(self, WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD):
        # constants
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.GENOME_LENGTH = GENOME_LENGTH
        self.LOCAL_WINDOW = LOCAL_WINDOW
        self.THRESHOLD = THRESHOLD
        self.COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]




# getting all the tag files for tag directory, filtering for necessary columns
# returns total number of tags
    def gather_data(self, data_dict, directory):
        tag_list = []
        tag_list = glob.glob(f"{directory}/*.tsv")
        cnt = 0
        for tag_file in tag_list:
            tag = pd.read_csv(tag_file, sep='\t', header=None)
            tag.columns = self.COLUMNS
            # removing unnecesary
            tag_filt = tag[['chromosome', 'position', 'read_len', 'strand']]
            data_dict[tag_file] = tag_filt
            cnt += len(tag_filt)
        return cnt


    # runs through all the tags for each sample
    # for each read in the sample updates counts based on overlap
    def get_counts(self, data, dictionary):
        for tag_file, tag_filt in data.items():
            for index, row in tag_filt.iterrows():
                self.overlap(dictionary, row['position'], row['read_len'], row['strand'])


    def get_dict_tags(self, df):    
        pos_list = list(df['position'])

        pos_dict = dict()

        for index in pos_list:
            pos_dict[index] = True
        return pos_dict  


    # overlap – increments all the values in the dictionary of windows 
    # tag – start index of tag, tag_len – length of tag, dictionary – read
    def overlap(self, dictionary, tag, tag_len, strand):
            center = int(tag + tag_len / 2)
            # the range of for loop is representing the starting positions of the window
            if strand: # checks for strand direction
                for x in range(center - self.WINDOW_LENGTH, center + 1):
                    # checking if the position exists; if yes – increment
                    if dictionary.get(x) == None:
                        dictionary[x] = 1
                    else:
                        dictionary[x] += 1
            else:
                # the range of for loop is representing the starting positions of the window
                # range is reversed
                for x in range(center + self.WINDOW_LENGTH, center + 1, -1):
                    # checking if the position exists; if yes – increment
                    if dictionary.get(x) == None:
                        dictionary[x] = 1
                    else:
                        dictionary[x] += 1