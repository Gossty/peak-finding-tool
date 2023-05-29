import glob 
import pandas as pd



class Formating:
    def __init__(self, WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD):
        # constants
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.GENOME_LENGTH = GENOME_LENGTH
        self.LOCAL_WINDOW = LOCAL_WINDOW
        self.THRESHOLD = THRESHOLD
        self.COLUMNS = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]
        self.COLUMNS_FILT = ['chromosome', 'position', 'read_len', 'strand']



# getting all the tag files for tag directory, filtering for necessary columns
# returns total number of tags
    def gather_data(self, directory):
        tag_list = []
        total_data = pd.DataFrame(columns=self.COLUMNS_FILT)
        tag_list = glob.glob(f"{directory}/*.tsv")

        for tag_file in tag_list:
            tag = pd.read_csv(tag_file, sep='\t', header=None)
            tag.columns = self.COLUMNS
            # removing unnecesary
            tag_filt = tag[self.COLUMNS_FILT]
            total_data = pd.concat([total_data, tag_filt])
        return total_data


    # runs through all the tags for each sample
    # for each read in the sample updates counts based on overlap
    def get_counts(self, dataframe, dictionary):
        for index, row in dataframe.iterrows():
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