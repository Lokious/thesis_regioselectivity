
import pandas as pd
import openpyxl
from sys import argv
def readrhlist():
    # rheauniprot = open("C:/Users/HP/Desktop/thesis/data/rhea2uniprot_sprot.tsv").readlines()
    # for line in rheauniprot:
    #     line_list = (line.strip()).split("\t")
    #     RHEA_ID = line_list[0]
    #     DIRECTION = line_list[1]
    #     uniprotid = line_list[3]
    #     if line_list[0].startswith("R"):
    #         rheauniprot_dataframe = pd.DataFrame(columns=[line_list[0], line_list[1]])
    #     else:
    #         new_row = {"RHEA_ID":RHEA_ID,"DIRECTION":DIRECTION}
    #         rheauniprot_dataframe = pd.concat([rheauniprot_dataframe, pd.DataFrame(new_row, index=[uniprotid])])
    rheauniprot_dataframe = pd.read_csv("C:/Users/HP/Desktop/thesis/data/rhea2uniprot_sprot.tsv", sep='\t', header=0, index_col=3)
    rheauniprot_dataframe = rheauniprot_dataframe.rename_axis(index='Entry')
    return rheauniprot_dataframe


def read_sequence(file):
    lines = pd.read_excel(file,index_col=0)
    return lines


def main():
    # readfile whihc contains the rhea id and related uniprotid
    rheauniprot_dataframe = readrhlist()
    print(rheauniprot_dataframe)
    #read id and sequences in dataframe
    id_seq_dataframe = read_sequence("C:/Users/HP/Desktop/thesis/data/id_tosequence.xlsx")
    print(id_seq_dataframe)
    #link the sequences and reaction participant put in one csv file
    fulldata = id_seq_dataframe.join(rheauniprot_dataframe)
    # print(fulldata)
    # groups_uniprotid = fulldata.groupby(["Entry"])
    # for group in groups_uniprotid:
    #     print(group)
    fulldata.to_csv("C:/Users/HP/Desktop/thesis/data/fulldata.csv")

if __name__ == "__main__":
    main()