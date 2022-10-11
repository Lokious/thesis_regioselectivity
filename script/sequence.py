#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
This code is to parse sequence data
"""
import os
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
import matplotlib.pyplot as plt

from os import path
import subprocess

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

    def remove_redundant_sequence(self,sequences_file="../autodata/rawdata/uniprot_ec2.1.1.fasta"):

        try:
            file=open("../autodata/sequences/uniprot_ec2.1.1_rm.fasta")
        except:
            print("removing redundant sequences from {}".format(sequences_file))
            if subprocess.check_output("bash runmmseqs.sh"):
                file = open("../autodata/sequences/uniprot_ec2.1.1_rep_seq.fasta")

    def get_active_site_dictionary_from_file(self,file="../autodata/align/align_seed_sequences_with_structure/3ROD_active_site_N.txt")->dict:
        """
        This function is to read active sits as dictionary

        :param file: path of tesxt file which saves the active sits information
        :return:active_sites_dictionary for example {11:"Y"}
        """
        try:
            acs_dataframe = pd.read_table(file, header=0, comment="#",
                                          index_col="auth")
            print("acs_dataframe")
            print(acs_dataframe.columns)
        except EOFError:
            print("please check file exist and in correct format")
        # for given index_col, either given as string name or column index.
        # If a sequence of int / str is given, a MultiIndex is used.

        active_sites_dictionary = {}
        print(acs_dataframe["AA"])
        for index in acs_dataframe.index:
            active_sites_dictionary[index]=acs_dataframe.loc[index,'AA']
            print(acs_dataframe.loc[index,'AA'])
        print(active_sites_dictionary)
        return active_sites_dictionary

    def get_AA_within_distance_from_structure_file(self,pdb_name="3rod.pdb",residue_dictionary:dict={},chain_id:str="A",group:str=""):
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
        file="../autodata/structure/{}".format(pdb_name)
        if residue_dictionary != {}:
            print("use input residue dictionary")
        else:
            try:
                print("get residue_dictionary from ../autodata/structure/{}_active_site_{}.txt ".format(pdb_name.split(".")[0],group))
                residue_dictionary=self.get_active_site_dictionary_from_file("../autodata/structure/{}_active_site_{}.txt".format(pdb_name.split(".")[0],group))
            except:
                # try:
                #     path=input("please input the path of active site file")
                #     residue_dictionary = self.get_active_site_dictionary_from_file(path)
                # except:
                    raise EOFError ("can not get residue_dictionary from get_active_site_dictionary_from_file function")
        # create parser
        parser = PDBParser()

        # read structure from file
        structure = parser.get_structure("pdb_structure",file)
        model = structure[0]
        chain = model[chain_id]
        residues = list(chain.get_residues())
        active_sites = []

        for id in residue_dictionary.keys():
            print(id)
            amino_acid = residue_dictionary[id]
            # check the id from the dictionary and id from pdb structure refers
            # to the same amino acid
            if protein_letters_1to3[amino_acid].upper() == chain[id].get_resname().upper():
                active_sites.append(chain[id])
            else:
                print("AA in active site file:{}, AA in pdb chainA: {}".format(protein_letters_1to3[amino_acid].upper(),chain[id].get_resname().upper()))
                print("Warning: The amino acid in active site file is not the same as that from pdb file")
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

        # structure_seq_length = 0
        for residue in residues:
            if residue.id[0] == ' ':
                #can not add to count the length because it does not start at 1
                # structure_seq_length += 1
                print(residue.id)
                print(residue.get_resname())
                structure_seq_length = residue.id[1]
            else:

                break
        print("length of the sequence: {}".format(structure_seq_length))
        return amino_acide_close_to_active_site, structure_seq_length

    def get_sites_from_alignment(self, fileformat="fasta", file="",
                                 active_site_dictionary: dict={},
                                 start_pos: int = 0,
                                 structure_chain="3rod.pdb_chainA_s001",group:str="",pdb_name="3rod.pdb") -> pd.DataFrame:
        """
        This is the function to get the site close to active site for aligned sequences

        :param fileformat:string, format of input alignment
        :param file:string, path and name of the aligned sequences
        :param start_pos: int, the starting position of the sequneces compare to the author annnotated index(start from 1)
        :return:

        """
        if active_site_dictionary=={}:
            print("Need active_site_dictionary, creating....")
            #split and remove the bit score and coverage
            group_1=group.split("_")[0]
            active_site_dictionary,sequence_length =self.get_AA_within_distance_from_structure_file(pdb_name=pdb_name,group=group_1)
            print(active_site_dictionary)


        #read the alignment and save to dataframe
        print("read alignment from: {}".format(file))
        align = AlignIO.read(file, fileformat)
        #print(align.get_alignment_length())
        align_array = np.array([list(rec) for rec in align])
        ids = list(rec.id for rec in align)
        align_pd = pd.DataFrame(data=align_array,index=ids)
        if fileformat !="fasta":
            align_pd.replace(".", "-", inplace=True)
        #remove sequences from seed
        for i in align_pd.index:
            if i != structure_chain:
                if "-" in i:
                    align_pd.drop(index=i,inplace=True)
                    ids.remove(i)

        print(align_pd)
        print(ids)
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
        print(align_pd)

        align_pd=align_pd.T
        align_pd.replace( np.nan,"-", inplace=True)
        print("length columns:{}".format(len(align_pd.columns)))
        print(len(align_pd.columns)+start_pos-1)
        print(sequence_length)
        #print('Q9KUA0' in align_pd.index)
        #  double check the sequences length are the same with the structure chain sequences
        if (len(align_pd.columns)+start_pos-1) != sequence_length:
            raise ValueError("The alignment length is not as same as the structure sequence")

        #print(align_pd)
        #create dataframe for saving the site close to active sites
        active_site_pd=pd.DataFrame(index=ids)
        #print(active_site_pd)

        print("only keep site close to active site")
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
        print("active_site pd:")
        print(active_site_pd)

        #drop sites with over 80% gap
        for column in active_site_pd.columns:
            value_count=active_site_pd[column].value_counts()
            if "-" in value_count:
                print("drop position which has over 50% gaps")
                if (value_count["-"]/len(active_site_pd.index)) > 0.5:
                    print("drop:{}".format(column))
                    active_site_pd.drop(columns=column,inplace=True)
        print(active_site_pd)

        count=0
        #construct MSA with amino acids close to active sites
        ##save sequences to fasta file
        file = open(
            "../autodata/protein_encoding/active_site/{}_seq_close_to_active_sites.fasta".format(
                group), "w")
        for seq_id in active_site_pd.index:
            seq= "".join(active_site_pd.loc[seq_id,:].values.tolist())
            #print(seq,seq_id)
            if "-" in seq:
                count +=1
            file.write(">{}\n".format(seq_id))
            file.write("{}\n".format(seq))
        print("count:{}".format(count))
        file.close()
        ##save pd.Dataframe to csv file
        active_site_pd.to_csv("../autodata/protein_encoding/active_site/{}_active_site_df.csv".format(group))
        return active_site_pd

    def amino_acid_properties(self, amino_acid: str, structure_aa: str) -> dict:
        """
        return list of amino acid properties generate by quantiprot

        :param
        amino_acid:three letter(all capital) or one letter code of amino acid
        structure_aa: amino acid from the samr position of structure sequence
        :return: properties_dictionary
        """

        structure_aa=structure_aa.upper()
        if len(amino_acid) == 3:
            amino_acid = protein_letters_3to1[amino_acid]
        else:
            amino_acid = amino_acid.upper()
        try:
            charge = get_aa2charge().mapping[amino_acid]
            volume = get_aa2volume().mapping[amino_acid]
            hydropathy = get_aa2hydropathy().mapping[amino_acid]
            hydrophobicity = get_aa2mj().mapping[amino_acid]
            if amino_acid == "-":
                similarity_score = bl.BLOSUM(62)[("*"+structure_aa)]
            else:
                print(amino_acid+structure_aa)
                similarity_score = bl.BLOSUM(62)[(amino_acid+structure_aa)]
            print(similarity_score)
            properties_dictionary = {"charge": charge, "volume": volume,
                                     "hydropathy": hydropathy,
                                     "hydrophobicity": hydrophobicity,
                                     "similarity_score": similarity_score}
        except:
            print("please check input {} in the amino acids list".format(amino_acid))
            properties_dictionary = {"charge": -1000, "volume": -1000,
                                     "hydropathy": -1000,
                                     "hydrophobicity": -1000,
                                     "similarity_score":-4}
        return properties_dictionary

    def drop_sequences(self,sequences_file="",length=None):
        """This function is to remove too long and too short sequences from fasta file

        :param sequences_file: string, name and path of fasta file
        :return: record_dict
        """

        record_dict= {rec.id: rec.seq for rec in SeqIO.parse(sequences_file, "fasta")}
        print("sequence number:{}".format(len(record_dict)))
        #save sequences length
        length_dictionary = {}
        for key in record_dict.keys():
            length = len(record_dict[key])
            #print(length)
            #print(record_dict[key])
            if length not in length_dictionary:
                length_dictionary[length] = 1
            else:
                length_dictionary[length] += 1
        if length ==None:
            sorted_length_dict = sorted(length_dictionary.items(), key=lambda x: x[1], reverse=True)
            #get the most frequent sequences length, use the first 10 to calculate the sverage length
            print(sorted_length_dict[:20])
            average=(sum([x[0] for x in sorted_length_dict[:20]])/20)
            print("average length of top 20 length{}".format(average))
        else:
            print("drop sequences whose length are 50 AA longer or shorter than structure sequences")
            average = length
        #remove too long and too short sequences
        remove_list = []
        for entry in record_dict.keys():
            length_seq=len(record_dict[entry])
            #if the distance between sequence length and average length larger than 100
            if abs(length_seq-average) > 100:
                print("remove length {}".format(length_seq))
                remove_list.append(entry)

        for entry in remove_list:
            del record_dict[entry]
        #print(len(record_dict))
        #write left sequences to file
        file=open(sequences_file,"w")
        for key in record_dict:
            # print(key)
            # print(record_dict[key])
            file.write(">{}\n".format(key))
            file.write("{}\n".format(record_dict[key]))
        print("Finish removing sequences!")
        return record_dict

    def remove_sequences_from_result_of_mmseqs(self,tab_file=".tab",seq_file=".fasta"):
        """
        remove high similarity sequences (0.9)get from mmseqs map

        :param tab_file:
        :param seq_file:
        :return:
        """
        tab_result= pd.read_csv(tab_file,sep="\t",header=None)
        tab_result=tab_result.iloc[:, :3]
        print(tab_result)
        record_dict= {rec.id: rec.seq for rec in SeqIO.parse(seq_file, "fasta")}
        #print(record_dict)
        print("sequence number:{}".format(len(record_dict)))

        #list of entry will be kept
        remove_list = []
        #list of entry will be removed
        keep_list = []
        ##NB: if entry1 and entry2 are both possiable to be kept keep entry1###
        for index in tab_result.index:
            entry1 = tab_result.loc[index,0]
            entry2 = tab_result.loc[index,1]
            if entry1 == entry2:
                continue
            #check if entry1 in remove list(similar to at least one entry in keep list)
            if entry1 not in remove_list:
                #check if entry1 in keep list(not similar to all entries in keep list)
                #one in keep, another should be add to remove list
                if entry1 in keep_list:
                    if entry2 not in keep_list:
                        remove_list.append(entry2)
                    else:
                        #print(entry1,entry2)
                        keep_list.remove(entry2)
                        remove_list.append(entry2)
                else:
                    if entry2 in keep_list:
                         remove_list.append(entry1)
                    else:
                        #if both of entry1 and 2 are not in keep list and entry1 not in remove list,
                        #add entry1 to keep and entry2 to remove
                        remove_list.append(entry2)
                        ###this is the one situation add entry1 to keep list
                        keep_list.append(entry1)
            else:
                #becaue if one entry is in remove list, it could not be add to keep,
                #so do not need to check
                if entry2 not in remove_list:
                    if entry2 in keep_list:
                        continue
                    else:
                        ###this is the one situation add entry2 to keep list
                        keep_list.append(entry2)
                else:
                    #both in remove, continue
                    continue


        keep_list = list(set(keep_list))
        remove_list=list(set(remove_list))

        for entry in remove_list:
            if entry in keep_list:
                print(entry)
                raise ValueError("OVERLAP BETWEEN KEEP AND REMOVE LIST!!")
            #del record_dict[entry]
        print(len(record_dict))
        #write left sequences to file
        file=open(seq_file,"w")
        for key in record_dict:
            # print(key)
            # print(record_dict[key])
            file.write(">{}\n".format(key))
            file.write("{}\n".format(record_dict[key]))
        print("Finish removing redundant sequences!")
        return record_dict
def main():

    #unittest.main()
    seq=Sequences()
    #seq.get_active_site_dictionary_from_file()
    #seq.get_AA_within_distance_from_structure_file()

    # seq.get_sites_from_alignment(
    #     file="../autodata/align/align_seed_sequences_with_structure/N_3rod_align_sequences",structure_chain="3rod.pdb_chainA_s001",pdb_name="3rod.pdb",
    #     start_pos=3,group="N")

    #seq.get_sites_from_alignment(file="../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/PF08241.15PF03602.18/PF08241.15PF03602.18_hmmalign_out_pdb_5WP4.aln",fileformat="clustal",start_pos=0, structure_chain="5WP4_1|Chain",group="PF08241.15PF03602.18",pdb_name="5wp4.pdb")
    seq.remove_sequences_from_result_of_mmseqs("../autodata/sequences/alnRes.tab",seq_file="../autodata/sequences/PF08241.15PF03602.18.fasta")
    #seq.group_seq_based_on_methylated_type(inputfile="../autodata/seq_smiles_all.csv",save_directory="../autodata/sequences")
    #seq.remove_duplicate()
    #seq.group_fg_based_on_methylated_type("../autodata/seq_smiles_all.csv",)

    #seq.drop_sequences("../autodata/sequences/PF08242.fasta")
if __name__ == "__main__":
    main()