import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        prog='peakFinding',
        description='Finds Peaks in the provided tag directory'
    )
    parser.add_argument('tag_directory', help="tag directory for analysis")
    parser.add_argument("control", help="control of peak finding")
    parser.add_argument('-style', help='<factor> for TFs; <histone> for histone modifications')

    args = parser.parse_args()
    print(args.control)
    print(args.tag_directory)
    gather_data(args.control, args.tag_directory)
    

def gather_data(control, tag_directory):
    chr1 = pd.read_csv(tag_directory, sep='\t')
    input1 = pd.read_csv(tag_directory, sep='\t')
    columns = ["blank", "chromosome", "position", "strand", "num_reads", "read_len"]
    chr1.columns = columns
    input1.columns = columns

    chr1_filt = chr1[['position', 'read_len']]
    input1_filt = input1[['position', 'read_len']]

    # check that it is actually sorted

    chr1_filt.sort_values(by=["position"], ascending=True, inplace=True)
    input1_filt.sort_values(by=["position"], ascending=True, inplace=True)

    
    print(input1_filt)


main()