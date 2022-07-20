#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
datafile should be put in the directory data/ ,this script should be under regioselectivity_prediction
includes Rhea-ec_2_1_1.tsv
"""

import pandas as pd
from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import pikachu
from script.molecular_class import molecular
from script import parse_data
import glob
def merge_uniprot_id_smile(rheauniprot_dataframe,seq_df):
    """
    combine sequence and substrate smile in onefile and save to csv file
    :param
        rheauniprot_dataframe:
        seq_df:
    :return:
    """
    df1 = parse_data.get_substrate_chebi("data/Rhea-ec_2_1_1.tsv")
    df1 = parse_data.get_smile(df1)
    df1.index = df1.index.map(lambda x: x.split(":")[1])
    df1["RHEA_ID"] = df1.index
    df1.set_index("RHEA_ID", inplace=True)
    uniprot_entry = parse_data.remove_duplicated_id(
        r"E:\Download\regioselectivity_prediction\data\hmm_out")
    fulldata = parse_data.rheaid_to_uniprot(uniprot_entry,
                                            rheauniprot_dataframe)
    fulldata.set_index("RHEA_ID", inplace=True)
    fulldata = fulldata.reset_index()
    df1 = df1.reset_index()
    complete_data = pd.merge(fulldata, df1)
    seq_df = seq_df.reset_index()
    complete_data = pd.merge(seq_df,complete_data)
    complete_data.to_csv("data/seq_smiles.csv")
    return complete_data
def main():
    # readfile whihc contains the rhea id and related uniprotid
    #run with command line
    # rh_file = argv[1]
    #run with pycharm
    rh_file = "data/rhea2uniprot_sprot.tsv"
    rheauniprot_dataframe = parse_data.readrhlist(rh_file)

    #read id and sequences in dataframe
    seq_file = "data/id_tosequence.xlsx" # run with pycharm
    #seq_file = argv[2]
    id_seq_dataframe = parse_data.read_sequence(seq_file)
    print(id_seq_dataframe)
    merge_uniprot_id_smile(rheauniprot_dataframe,id_seq_dataframe)




    #print(id_seq_dataframe)
    #link the sequences and reaction participant put in one csv file
    #fulldata = id_seq_dataframe.join(rheauniprot_dataframe)
    # print(fulldata)
    # groups_uniprotid = fulldata.groupby(["Entry"])
    # for group in groups_uniprotid:
    #     print(group)
    #save to file
    # file_name = input("please input file name and directory")
    # fulldata.to_csv(file_name)
    #convert rheaid to uniprot entry, as there are no direct link between Rhea id to sequences
    #read rhea id file
    #rhea = open("data/Rhea-ec_2.1.1.list").readlines()
    #uniprot_list = rheaid_touniprot(rhea, rheauniprot_dataframe)
    #print(uniprot_list)
    # file = open("uniprot_list.txt","w")
    # for index,entry in enumerate(uniprot_list):
    #     if isinstance(entry, str):
    #         file.write("{},{}\n".format(entry,rhea[index]))
    #     else:
    #         entry.to_csv("uniprot_list.txt", mode='a',header=False)
    # print(uniprot_list)


    #if the nmbering from rdkit is following the carbon numbering rules?
    # mol = mol_with_atom_index('CC1=C[N]C2=C1[C@H](C)C=CC2=O.C(C5=CC4=C3CCN(C3=C(C(=C4[N]5)OC)O)C(=O)N)(=O)N6C7=C(CC6)C8=C(C(=C7O)OC)[N]C(=C8)C(=O)[N]9C%10=C(C=C9)C%11=C(C(=C%10)O)N=CC%11C')
    # Draw.ShowMol(mol,size=(600,600))
    # smiles_with_atom_mappings = Chem.MolToSmiles(mol)



    #remove_duplicated_id("data/hmm_out/top_ten_hits_exclude_nomethylrelated/hmmscantbout_top_6_hit.pfam")

    # uniprot_entry = remove_duplicated_id(r"E:\Download\regioselectivity_prediction\data\hmm_out")
    # rheaid_to_uniprot(uniprot_entry, rheauniprot_dataframe)
    #pikachu.general.draw_smiles(
    #    "CCC(=O)[C@@H]1NC(=O)[C@H](CC(N)=O)NC(=O)C(NC(=O)[C@@H](N)CC(C)C)[C@H](O)c3ccc(Oc2cc1cc(CC)c2O)c(Cl)c3")
if __name__ == "__main__":
    main()
