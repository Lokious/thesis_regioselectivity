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
from script.molecular_class import molecular, reaction
from script import parse_data

import pandas as pd
from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
import pikachu
from rdkit.Chem import Draw, AllChem
import glob
from IPython.display import SVG, display, Image
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from PIL import Image

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
    complete_data = parse_data.merge_reaction(complete_data)
    complete_data.to_csv("data/seq_smiles.csv")
    return complete_data


def drop_useless_column(dataframe):
    dataframe = dataframe.loc[:,["Sequence","sub_mols","pro_mols","sub_smiles","pro_smiles"]]
    #print(dataframe)


def keep_longest_smile(dataframe):

    dataframe["main_sub"] = pd.DataFrame(
        len(dataframe.index) * [0])
    dataframe["main_pro"] = pd.DataFrame(
        len(dataframe.index) * [0])
    for index in dataframe.index:
        main_sub = max((dataframe.loc[index,"sub_smiles"]).split("."), key=len)
        dataframe.loc[index, "main_sub"] = main_sub
        main_pro = max((dataframe.loc[index,"pro_smiles"]).split("."), key=len)
        dataframe.loc[index, "main_pro"] = main_pro
    return dataframe


def return_reactions(dataframe):
    for index in dataframe.index:
        rxn = dataframe.loc[index,"rxn"]
        sub = dataframe.loc[index, "main_sub"]
        pro = dataframe.loc[index, "main_pro"]
        reaction1 = reaction(substrates=sub, products=pro)
        #print(dataframe.loc[index,"RHEA_ID"])
        r1 = reaction1.get_reaction_sites(rxn_object=rxn)

        r2 = reaction1.get_reactant_atom()
        break



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
    seq_smiles = merge_uniprot_id_smile(rheauniprot_dataframe,id_seq_dataframe)
    data_frame = keep_longest_smile(seq_smiles)

    return_reactions(seq_smiles)
    drop_useless_column(seq_smiles)



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
    # mol = molecular_class.mol_with_atom_index('CC1=C[N]C2=C1[C@H](C)C=CC2=O.C(C5=CC4=C3CCN(C3=C(C(=C4[N]5)OC)O)C(=O)N)(=O)N6C7=C(CC6)C8=C(C(=C7O)OC)[N]C(=C8)C(=O)[N]9C%10=C(C=C9)C%11=C(C(=C%10)O)N=CC%11C')
    # Draw.ShowMol(mol,size=(600,600))
    # smiles_with_atom_mappings = Chem.MolToSmiles(mol)



    #remove_duplicated_id("data/hmm_out/top_ten_hits_exclude_nomethylrelated/hmmscantbout_top_6_hit.pfam")

    # uniprot_entry = remove_duplicated_id(r"E:\Download\regioselectivity_prediction\data\hmm_out")
    # rheaid_to_uniprot(uniprot_entry, rheauniprot_dataframe)
    #pikachu.general.draw_smiles(
    #    "cc")
    # rxn = AllChem.ReactionFromSmarts(r"[C:1]-[C:2]>>[C:1].[C:2]")
    # img = Draw.ReactionToImage(rxn,returnPNG=True)
    # #?
    # with open("test.png",'wb') as handle:
    #     handle.write(img)


if __name__ == "__main__":
    main()
