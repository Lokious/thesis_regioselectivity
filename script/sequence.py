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

from Bio import AlignIO, SeqIO
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_1to3

from os import path

class Sequences():

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

    def get_sites_from_alignment(self,format="fasta",file=""):
        """


        :param format:
        :param file:
        :return:


        """
        align = AlignIO.read(file, "fasta")

    def get_AA_within_distance_from_structure_file(self,file="../autodata/align/align_seed_sequences_with_structure/3rod.pdb",residue_dictionary:dict={11:"Y"}):
        """
        This function is to get the id of residues close to active sites

        :param file: string, path and name of pdb file
        :param residue_dictionary: {id of active site:only letter code of AA}
        :return:
        amino_acide_close_to_active_site: dictionary, id of active site is key,
        value is a tuple which contains the id, name and distance of close AA
        For example:
        {11: [(8, 'LYS', 5.5259695), (10, 'THR', 3.8049822), (11, 'TYR', 0.0)}


        Some related problem: https://www.biostars.org/p/401105/
        """

        # create parser
        parser = PDBParser()

        # read structure from file

        structure = parser.get_structure("pdb_structure",file)
        model = structure[0]
        chain = model['A']
        residues=list(chain.get_residues())
        active_sites=[]
        for id in residue_dictionary.keys():
            amino_acid= residue_dictionary[id]
            # check the id from the dictionary and id from pdb structure refers
            # to the same amino acid
            if protein_letters_1to3[amino_acid].upper() == chain[id].get_resname().upper():
                active_sites.append(chain[id])

        amino_acide_close_to_active_site={}
        for id1,residue1 in enumerate(active_sites):
            #remove residue from the list of residues, so it only count the distance between two residue once


            #seems residue.get_id()[1] will return the id from author, which is
            #use in the paper

            for id2,residue2 in enumerate(residues):
                # compute distance between alpha C atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom
                    continue
                #save AA to dictionary whihc are close to the active site
                if distance <= 5:
                    if residue1.get_id()[1] not in amino_acide_close_to_active_site:
                        amino_acide_close_to_active_site[residue1.get_id()[1]]=[(residue2.get_id()[1],residue2.get_resname(),distance)]
                    else:
                        amino_acide_close_to_active_site[residue1.get_id()[1]].append((residue2.get_id()[1],residue2.get_resname(),distance))

        print(amino_acide_close_to_active_site)
        return amino_acide_close_to_active_site

def main():
    #unittest.main()

    seq=Sequences()
    seq.get_AA_within_distance_from_structure_file()

    #seq.group_seq_based_on_methylated_type(inputfile="../autodata/seq_smiles_all.csv",save_directory="../autodata/sequences")
    #seq.remove_duplicate()
    # seq.group_fg_based_on_methylated_type("../autodata/seq_smiles_all.csv",)

if __name__ == "__main__":
    main()