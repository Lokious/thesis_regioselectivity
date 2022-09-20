#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
This code is to parse sequence data
"""
import unittest
import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_1to3,protein_letters_3to1
from quantiprot.metrics.aaindex import get_aa2volume, get_aa2hydropathy, get_aa2charge,get_aa2mj
import blosum as bl


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


    def get_active_site_dictionary_from_file(self,file="../autodata/align/align_seed_sequences_with_structure/3ROD_active_site.txt")->dict:
        """
        This function is to read active sits as dictionary

        :param file: path of tesxt file which saves the active sits information
        :return:active_sites_dictionary for example {11:"Y"}
        """
        try:
            acs_dataframe = pd.read_table(file, header=0, comment="#",
                                          index_col="auth")
        except EOFError:
            print("please check file exist and in correct format")
        # for given index_col, either given as string name or column index.
        # If a sequence of int / str is given, a MultiIndex is used.

        active_sites_dictionary={}
        for index in acs_dataframe.index:
            active_sites_dictionary[index]=acs_dataframe.loc[index,'AA']

        print(active_sites_dictionary)
        return active_sites_dictionary

    def get_AA_within_distance_from_structure_file(self,file="../autodata/align/align_seed_sequences_with_structure/3rod.pdb",residue_dictionary:dict={11:"Y"},chain_id:str="A"):
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
        try:
            residue_dictionary=self.get_active_site_dictionary_from_file()
            print(residue_dictionary)
        except:
            print("can not get residue_dictionary from get_active_site_dictionary_from_file function")
        # create parser
        parser = PDBParser()

        # read structure from file
        structure = parser.get_structure("pdb_structure",file)
        model = structure[0]
        chain = model[chain_id]
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
            for id2,residue2 in enumerate(residues):

                # compute distance between alpha C atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom
                    continue
                #save AA to dictionary whihc are close to the active site
                if distance <= 5:
                    # seems residue.get_id()[1] will return the id from author, which is
                    # use in the paper
                    if residue1.get_id()[1] not in amino_acide_close_to_active_site:
                        amino_acide_close_to_active_site[residue1.get_id()[1]]=[(residue2.get_id()[1],residue2.get_resname(),distance)]
                    else:
                        amino_acide_close_to_active_site[residue1.get_id()[1]].append((residue2.get_id()[1],residue2.get_resname(),distance))

        structure_seq_length = 0
        for residue in residues:
            if residue.id[0] == ' ':
                structure_seq_length += 1
                print(residue.id)
                print(residue.get_resname())
        print(structure_seq_length)
        return amino_acide_close_to_active_site, structure_seq_length

    def get_sites_from_alignment(self, fileformat="fasta", file="",
                                 active_site_dictionary: dict={},
                                 start_pos: int = 0,
                                 structure_chain="3rod.pdb_chainA_s001") -> pd.DataFrame:
        """
        This is the function to get the site close to active site for aligned sequences

        :param fileformat:string, format of input alignment
        :param file:string, path and name of the aligned sequences
        :param start_pos: int, the starting position of the sequneces compare to the author annnotated index
        :return:

        """
        if active_site_dictionary=={}:
            print("Need active_site_dictionary, creating....")
            active_site_dictionary,sequence_length =self.get_AA_within_distance_from_structure_file()
            print(active_site_dictionary)


        #read the alignment and save to dataframe
        align = AlignIO.read(file, fileformat)
        print(align.get_alignment_length())
        align_array = np.array([list(rec) for rec in align])
        ids = list(rec.id for rec in align)
        align_pd = pd.DataFrame(data=align_array,index=ids)
        print('Q9KUA0' in list(align_pd.index))

        #drop the column which the guided structure is a gap, and drop 'X' amino acid
        ## change - and X
        align_pd.replace("X", "-", inplace=True)
        align_pd.replace("-", np.nan, inplace=True)
        align_pd = align_pd.T
        align_pd.dropna(subset=[structure_chain],inplace=True)
        align_pd.reset_index(drop=True, inplace=True)

        #rest column names to fit the index from author
        index_list=list(align_pd.index)
        align_pd.index=[int(x)+start_pos for x in index_list]
        print(len(align_pd.columns))
        align_pd=align_pd.T
        align_pd.replace( np.nan,"-", inplace=True)
        print('Q9KUA0' in align_pd.index)
        #  double check the sequences length are the same with the structure chain sequences
        if len(align_pd.columns) != sequence_length:
            raise ValueError("The alignment length is not as same as the structure sequence")

        #print(align_pd)
        #create dataframe for saving the site close to active sites
        active_site_pd=pd.DataFrame(index=ids)
        #print(active_site_pd)
        for key_acs in active_site_dictionary.keys():
            for item_tuple in active_site_dictionary[key_acs]:
                #activesite_closeAA
                columnname=item_tuple[0]
                active_site_pd[columnname]=["NA"]*len(ids)

        for id in ids:
            for column in active_site_pd.columns:
                #print(column)
                aa = align_pd.loc[id,column]
                active_site_pd.loc[id,column]=aa
        print(active_site_pd)

        #construct MSA with amino acids close to active sites
        ##save sequences to fasta file
        for seq_id in active_site_pd.index:
            file = open("../autodata/align/align_seed_sequences_with_structure/N_seq_close_to_active_sites.fasta", "a")
            seq= "".join(active_site_pd.loc[seq_id,:].values.tolist())
            print(seq,seq_id)
            file.write(">{}\n".format(seq_id))
            file.write("{}\n".format(seq))

        ##save pd.Dataframe to csv file
        active_site_pd.to_csv("../autodata/align/align_seed_sequences_with_structure/N_active_site_df.csv")
        return active_site_pd

    def amino_acid_properties(self, amino_acid: str, structure_aa: str) -> dict:
        """
        return list of amino acid properties generate by quantiprot

        :param
        amino_acid:three letter(all capital) or one letter code of amino acid
        structure_aa: amino acid from the samr position of structure sequence
        :return: properties_dictionary
        """
        if len(amino_acid) == 3:
            amino_acid = protein_letters_3to1[amino_acid]
        try:
            charge = get_aa2charge().mapping[amino_acid]
            volume = get_aa2volume().mapping[amino_acid]
            hydropathy = get_aa2hydropathy().mapping[amino_acid]
            hydrophobicity = get_aa2mj().mapping[amino_acid]
            similarity_score = bl.BLOSUM(62)[(amino_acid+structure_aa)]
            print(similarity_score)
            properties_dictionary = {"charge": charge, "volume": volume,
                                     "hydropathy": hydropathy,
                                     "hydrophobicity": hydrophobicity,
                                     "similarity_score": similarity_score}
        except:
            print("please check input {} in the amino acids list".format(amino_acid))
            properties_dictionary = {"charge": 0, "volume": 0,
                                     "hydropathy": 0,
                                     "hydrophobicity": 0,
                                     "similarity_score":0}
        return properties_dictionary

def main():

    #unittest.main()
    seq=Sequences()
    #seq.get_active_site_dictionary_from_file()
    #seq.get_AA_within_distance_from_structure_file()

    #seq.get_sites_from_alignment(file="../autodata/align/align_seed_sequences_with_structure/N_3rod_align_sequences",start_pos=3)
    seq.amino_acid_properties("Y", "Y")
    #seq.group_seq_based_on_methylated_type(inputfile="../autodata/seq_smiles_all.csv",save_directory="../autodata/sequences")
    #seq.remove_duplicate()
    #seq.group_fg_based_on_methylated_type("../autodata/seq_smiles_all.csv",)


if __name__ == "__main__":
    main()