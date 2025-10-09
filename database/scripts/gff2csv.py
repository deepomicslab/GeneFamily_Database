import os
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from pathlib import PosixPath

def arg_parser():
    parser = ArgumentParser(
            formatter_class=ArgumentDefaultsHelpFormatter,
            conflict_handler="resolve", prog="gff2csv")
    parser.add_argument(
            "-i", "--hmmin", required=True, dest="gffinput",
            metavar="STR", type=str,
            help="input gene_exon_info.csv")
    parser.add_argument(
            "-o", "--output", required=True, dest="out",
            metavar="STR", type=str,
            help="output csv files"
    )
    return parser.parse_args()


def gff2csv(infile, outfile):
    change_list = []
    with open(infile, 'rt') as f1:
        for eachline in f1:
            array = eachline.split('\t')
            array[-1] = array[-1].replace('\n','')
            newarray = ''
            for i, element in enumerate(array):
                if i != len(array):
                    newarray = newarray + element + ','
                else:
                    newarray = newarray + element
            change_list.append(newarray + '\n')

    with open(outfile, 'wt') as f2:
        for i in change_list:
            f2.write(i)


def search_fa(p,outpath):
    for infile in p.iterdir():
        if infile.suffix == '.gff':
            outfile = Path(outpath) / (infile.stem + ".csv")
            gff2csv(infile, outfile)


def main():
    args = arg_parser()
    p = Path(args.gffinput)
    print("p",p)
    search_fa(p,args.out)

if __name__ == "__main__":
    main()
