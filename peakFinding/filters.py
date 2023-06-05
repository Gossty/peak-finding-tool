from scipy.stats import poisson
import numpy as np

# Class containing all filters used in the peakFinding tool
class Filters():

    # constructor for global values
    def __init__(self, WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                sample_length, control_length, FOLD_VALUE):
        # constants
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.GENOME_LENGTH = GENOME_LENGTH
        self.LOCAL_WINDOW = LOCAL_WINDOW
        self.THRESHOLD = THRESHOLD
        self.sample_length = sample_length
        self.control_length = control_length
        self.FOLD_VALUE = FOLD_VALUE


    # filtering based of fold change.
    # tf_bound keeps track of positions where TFs bound based on whether fold change < 4
    def fc_filt(self, sample_counts, control_counts, peaks_arr):
        tf_bound = []
        for index in peaks_arr:
            # check if control has this window
            if control_counts.get(index) == None:
                continue
            # normalizing and getting fold value
            fold_value = (sample_counts[index] / self.sample_length) / (
                          control_counts[index] / self.control_length)

            # checking fold value past threshold and adding to array
            if  fold_value >= self.FOLD_VALUE:
                tf_bound.append(index)
        tf_bound.sort()
        return tf_bound


    # goes through tf_bound and filters to position with max tag count in each window
    # window is defined by WINDOW_LENGTH * 2
    # returns an array of filtered starting positions
    def max_count_filt(self, tf_bound, sample_counts):
        max_filt = []
        if len(tf_bound) == 0:
            return max_filt
        check = tf_bound[0]
        # iterating through starting position of tf_bound
        for index in tf_bound:

            # continue iterating until the check
            if index < check:
                continue
            
            check = index       # check keeps track of the current window
            start = index       # starting value of for loop below

            position = start    # position represents max peak in a window
            max_dict = dict()   # helper dictionary for finding max peak in a window

            # iterating through the current window
            for i in range(start, start + (self.WINDOW_LENGTH * 2) + 1):
                # checking if the starting position of a window exists
                if sample_counts.get(i) == None:
                    continue   
                max_dict[i] = sample_counts.get(i)
            
            position = max(max_dict)

            # figure out the biggest peak
            max_filt.append(position)
            # we update check based on whether we passed the window
            check += (self.WINDOW_LENGTH * 2)

        return max_filt


    # dictionary for helper method
    def local_filt(self, sample_counts, tf_bound, dictionary):
        check = 0           
        start_check = 0    
        loc_density = 0     
        tf_bound_local_filt = []
        for index in tf_bound:

            if check <= 0:
                # local density defined by LOCAL_WINDOW
                loc_density = self.local_density(index, dictionary)
                check = self.LOCAL_WINDOW   # check keeps track of the current window
                prev_index = index          # previous index in the tf_bound

            # current desnity defined by WINDOW_LENGTH
            curr_density = sample_counts[index] / self.WINDOW_LENGTH 
            fold_value = curr_density / loc_density

            # checking fold value past threshold and adding to array
            if  fold_value >= self.FOLD_VALUE:
                tf_bound_local_filt.append(index)
            
            # check decrements based on the delta between prev and curr indexes
            check -= (index - prev_index)
            prev_index = index

        return tf_bound_local_filt

        # TODO don't count for tags that are inside of bigger peak
        
        # calculates density based on LOCAL_WINDOW
    def local_density(self, index, dictionary):
        density = 1

        # iterates through the current window
        for i in range(index, index + self.LOCAL_WINDOW + 1):
            if dictionary.get(i) == None:   # checking if a tag exists at given position
                continue
            density += 1
        density /= self.LOCAL_WINDOW
        return density


    # calculates p-value based on expected value for a given window and filters
    # by the given threshold
    def poisson_filt(self, tf_bound, sample_counts, tags_in_peaks=None):
        peaks = []
        #lambda for poisson/expected value
        exp = (self.WINDOW_LENGTH * self.control_length) / self.GENOME_LENGTH

        for index in tf_bound:
            # checking if the starting position of a window exists
            if sample_counts.get(index) == None:
                continue
            # calculating p_value based on cdf
            p_value =  1 - poisson.cdf(sample_counts[index], exp)

            # checking p-value past threshold and adding to array
            if p_value < self.THRESHOLD:
                peaks.append(index)
                if tags_in_peaks is not None:
                    tags_in_peaks += sample_counts[index]
        if tags_in_peaks is not None:
            return peaks, tags_in_peaks
        
        return peaks

