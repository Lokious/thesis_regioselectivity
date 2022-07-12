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

def rheaid_touniprot(rhea,data_frame):

    entrys = []
    for id in rhea:

        number = (id.strip()).split(":")[-1]

        uniprot_entry = data_frame[data_frame['RHEA_ID'] == int(number)]
        if uniprot_entry.empty:
            entrys.append("NA")
        else:

            entrys.append(uniprot_entry.iloc[:,0])
    return entrys
def get_substrate_inchkey():
    """
    This function is used to split the substrate and product and save in dataframe
    :parameter

    :return: datframe which use rhea as index
    """
    #read file includes rhea id and chEBI id and name
    rhea_dataframe = pd.read_table("data/Rhea-ec_2_1_1.tsv", header=0, index_col=0,
                                          sep='\t')
    rhea_dataframe["substrate"] = 0
    for id in list(rhea_dataframe.index):
        equation = rhea_dataframe.loc[id,"Equation"]
        subtrates = equation.split("=")[0]
        subtrate_list = subtrates.split("+")
        name = (rhea_dataframe.loc[id,"ChEBI name"]).split(";")
        chebiid = (rhea_dataframe.loc[id,"ChEBI identifier"]).split(";")
        participentid = (rhea_dataframe.loc[id,"Participant identifier"]).split(";")
        subtrate_dict_1 = dict(zip(name,chebiid))
        subtrate_dict_2 = dict(zip(name, participentid))
        print(subtrate_dict_1)
        sub = ""
        for item in subtrate_list:
            if item.strip(" ") in subtrate_dict_1:
                sub += subtrate_dict_1[item.strip(" ")]+";"
            else:
                #error may occur because some of the substrate is not exact
                sub += "error;"

        rhea_dataframe.loc[id,"substrate"] = sub
    print(rhea_dataframe)

def main():
    # readfile whihc contains the rhea id and related uniprotid
    #run with command line
    # rh_file = argv[1]
    #run with pycharm
    rh_file = "data/rhea2uniprot_sprot.tsv"
    rheauniprot_dataframe = readrhlist(rh_file)

    #print(rheauniprot_dataframe)
    #read id and sequences in dataframe
    # run with pycharm
    seq_file = "data/id_tosequence.xlsx"
    #seq_file = argv[2]
    id_seq_dataframe = read_sequence(seq_file)
    #print(id_seq_dataframe)
    #link the sequences and reaction participant put in one csv file
    fulldata = id_seq_dataframe.join(rheauniprot_dataframe)
    # print(fulldata)
    # groups_uniprotid = fulldata.groupby(["Entry"])
    # for group in groups_uniprotid:
    #     print(group)
    #save to file
    # file_name = input("please input file name and directory")
    # fulldata.to_csv(file_name)
    #convert rheaid to uniprot entry, as the directly find enzyme from Rhea seems doesn't work
    #read rhea id file
    # rhea = open("data/Rhea-ec_2.1.1.list").readlines()
    # uniprot_list = rheaid_touniprot(rhea, rheauniprot_dataframe)
    # file = open("uniprot_list.txt","w")
    # for index,entry in enumerate(uniprot_list):
    #     if isinstance(entry, str):
    #         file.write("{},{}\n".format(entry,rhea[index]))
    #     else:
    #         entry.to_csv("uniprot_list.txt", mode='a',header=False)
    # print(uniprot_list)

    get_substrate_inchkey()
if __name__ == "__main__":
    main()
