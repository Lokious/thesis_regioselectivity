#!/usr/bin/env python3
import pandas as pd
from sys import argv
def readrhlist(file_name):

    rheauniprot_dataframe = pd.read_table(file_name, header=0, index_col=3,sep='\t', encoding="ISO-8859-1")
    rheauniprot_dataframe = rheauniprot_dataframe.rename_axis(index='Entry')
    return rheauniprot_dataframe


def read_sequence(file):
    lines = pd.read_excel(file,index_col=0)
    return lines


def main():
    # readfile whihc contains the rhea id and related uniprotid
    rh_file = argv[1]
    rheauniprot_dataframe = readrhlist(rh_file)
    print(rheauniprot_dataframe)
    #read id and sequences in dataframe
    seq_file = argv[2]
    id_seq_dataframe = read_sequence(seq_file)
    print(id_seq_dataframe)
    #link the sequences and reaction participant put in one csv file
    fulldata = id_seq_dataframe.join(rheauniprot_dataframe)
    # print(fulldata)
    # groups_uniprotid = fulldata.groupby(["Entry"])
    # for group in groups_uniprotid:
    #     print(group)
    #save to file
    file_name = input("please input file name and directory")
    fulldata.to_csv(file_name)

if __name__ == "__main__":
    main()
