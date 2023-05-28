from scipy.stats import poisson


class Filters():


    def __init__(self, WINDOW_LENGTH, GENOME_LENGTH, LOCAL_WINDOW, THRESHOLD,
                sample_length, control_length):
        # constants
        self.WINDOW_LENGTH = WINDOW_LENGTH
        self.GENOME_LENGTH = GENOME_LENGTH
        self.LOCAL_WINDOW = LOCAL_WINDOW
        self.THRESHOLD = THRESHOLD
        self.sample_length = sample_length
        self.control_length = control_length


    # filtering based of log fold change.
    # tf_bound keeps track of positions where TFs bound based on whether log2 fold change < 4
    def fc_filt(self, sample_counts, control_counts):
        tf_bound = []
        for index in sample_counts.keys():
            # check if control has this window
            if control_counts.get(index) == None:
                continue
            fold_value = (sample_counts[index] / self.sample_length) / (control_counts[index] / self.control_length)
            if  fold_value >= 4:
                tf_bound.append(index)
                # tf_bound[index] = log_value
        tf_bound.sort()
        return tf_bound


    # goes through tf_bound and filters to position with max tag count in each window
    # window is defined by WINDOW_LENGTH * 2
    # returns an array of filtered starting positions
    def max_count_filt(self, tf_bound, sample_counts):
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
            for i in range(start, start + (self.WINDOW_LENGTH * 2) + 1):
                if sample_counts.get(i) == None:
                    continue   
                max_dict[i] = sample_counts.get(i)
            
            position = max(max_dict)
            # figure out the biggest peak
            max_filt.append(position)
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
                loc_density = self.local_density(index, dictionary)
                check = self.LOCAL_WINDOW
                start_check = index

            curr_density = sample_counts[index] / self.WINDOW_LENGTH
            fold_value = curr_density / loc_density
            if  fold_value >= 4:
                tf_bound_local_filt.append(index)
            
            check -= (index - start_check)
            start_check = index
        print("length of local_filt: ", len(tf_bound_local_filt))
        return tf_bound_local_filt


        # calculates a tag density for a 10,000 bp window
    def local_density(self, index, dictionary):
        density = 1

        for i in range(index, index + self.LOCAL_WINDOW + 1):
            if dictionary.get(i) == None:
                continue
            density += 1
        density /= self.LOCAL_WINDOW
        return density


    # calculates p-value based on poisson distribution for a given window and filters
    # by the given threshold
    def poisson_filt(self, tf_bound, sample_counts):
        peaks = []
        #lambda for poisson/expected value
        exp = (self.WINDOW_LENGTH * self.control_length) / self.GENOME_LENGTH
        # expected value – shows mean tags for input
        for index in tf_bound:
            if sample_counts.get(index) == None:
                continue
            p_value =  1 - poisson.cdf(sample_counts[index], exp)
            if p_value < self.THRESHOLD:
                peaks.append(index)

        print("total number of peaks in poisson for control: ", len(peaks))

        return peaks