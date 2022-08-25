#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
This code is to parse sequence data
"""
import sys
import unittest
import pandas as pd
from Bio import SeqIO
from os import path
class sequences():

    def __init__(self, msa_model=None):
        self.model = msa_model

    def group_seq_based_on_methylated_type(self,inputfile="../data/seq_smiles_all_MANUAL.csv",save_directory="../data/sequences"):
        """

        :param inputfile:
        :param save_directory:
        :return:
        """
        whole_data = pd.read_csv(inputfile,header=0,index_col=0)

        seq_type_df = pd.DataFrame()
        seq_type_df["Sequence"] = whole_data["Sequence"]
        seq_type_df["methyl_type"] = whole_data["methyl_type"]
        seq_type_df["Entry"]= whole_data["Entry"]
        seprate_dataset = seq_type_df.groupby(by=["methyl_type"])
        files=[]
        for group in seprate_dataset.groups:
            sub_df = seprate_dataset.get_group(group)
            group = sub_df["methyl_type"].unique()[0]
            print(group)
            files.append(group)
            if path.exists("{}/{}.fasta".format(save_directory,group)):
                print("{}/{}.fasta is already exit".format(save_directory,group))
                continue
            else:
                for index in sub_df.index:
                    file = open("{}/{}.fasta".format(save_directory,group),"a")
                    file.write(">{}\n".format(sub_df.loc[index,"Entry"]))
                    file.write("{}\n".format(sub_df.loc[index,"Sequence"]))
                else:
                    #close the file when fiinihsed writing
                    file.close()
    def group_fg_based_on_methylated_type(self,inputfile,numbit:int=2048,bond:int=2):
        """

        :param inputfile:
        :param numbit:
        :param bond:
        :return:
        """
        input_df = pd.read_csv(inputfile,header=0,index_col=0)
        seprate_dataset = input_df.groupby(by=["methyl_type"])
        for group in seprate_dataset.groups:
            sub_df = seprate_dataset.get_group(group)
            group = sub_df["methyl_type"].unique()
            print(group)
            sub_df.reset_index(drop=True, inplace=True)
            sub_df.to_csv(
                "../autodata/group/input_dataframe_dropatoms_{}_withtype_bond{}_{}.csv".format(numbit,bond,group))

    def remove_duplicate(self,files):
        for file in files:
            with open('../autodata/sequences/{}_rm.fasta'.format(file), 'a') as outFile:
                record_ids = list()
                for record in SeqIO.parse("../autodata/sequences/{}.fasta".format(file), 'fasta'):
                    if record.id not in record_ids:
                        record_ids.append(record.id)
                        SeqIO.write(record, outFile, 'fasta')
        else:
            print("finished removing duplicate")
    def rename_sequences(self):

        from Bio.SeqRecord import SeqRecord
        for file in ["O", "N", "C", "O_N", "S", "Co", "Se"]:
            with open('../autodata/sequences/all_rename_type.fasta'.format(file),
                      'a') as outFile:
                for record in SeqIO.parse("../autodata/sequences/{}_rm.fasta".format(file), 'fasta'):
                    record_rename=SeqRecord(
                    seq=record.seq,
                    id=(file),
                    name=record.name)
                    SeqIO.write(record_rename, outFile, 'fasta')
                    # print(record)
                    # print(record_rename)


def main():
    #unittest.main()

    seq=sequences()

    #seq.group_seq_based_on_methylated_type(inputfile="../autodata/seq_smiles_all.csv",save_directory="../autodata/sequences")
    #seq.remove_duplicate()
    # seq.group_fg_based_on_methylated_type("../autodata/seq_smiles_all.csv",)

if __name__ == "__main__":
    main()