import argparse

def main():
    parser = argparse.ArgumentParser(
        prog='Peak Finding Tool',
        description='Finds Peaks in the provided tag directory'
    )
    parser.add_argument('tag directory', help="tag directory for analysis")
    parser.add_argument("control", help="control of peak finding")
    parser.add_argument('-style', help='<factor> for TFs; <histone> for histone modifications')

    args = parser.parse_args()
    print(args.accumulate(args.integers))

main()