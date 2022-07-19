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
import glob
def target_sequences(file):
    """
    This function is used to return Uniprot entry of target sequences from hmmscan
    :param file:string, file name
    :return: hmmscan_df: dataframe, include all the Uniprot entry from hmmscan output
    """
    name = ['target_name', 'accession_domain','query_name', 'accession',
            'E-value1',  'score1',  'bias1',   'E-value',  'score',  'bias',
            'exp', 'reg', 'clu',  'ov', 'env', 'dom','rep','inc', 'description_of_target']
    hmmscan_df = (pd.read_table(file, header=None, comment='#',names=name ,sep= '\s+', skip_blank_lines=True)).dropna()
    hmmscan_df = hmmscan_df.iloc[:,2]
    for i in hmmscan_df.index:
        hmmscan_df.iloc[i] = (hmmscan_df.iloc[i]).split("|")[1]
    hmmscan_df.columns = ["index","id"]
    return hmmscan_df
def remove_duplicated_id(directory):

    file_name = glob.iglob(directory+'/*.tsv')
    datafrmae_list = []
    for i in file_name:
        df = target_sequences(i)
        print(df.columns)
        #????
        print(df.index)
        #print(df.iloc[:,1])
        datafrmae_list.append(df['id'])
        print(datafrmae_list)
    
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


def get_substrate_chebi(file):
    """
    This function is used to split the substrate and product and save in dataframe
    :parameter:None

    :return:
    rhea_dataframe: datframe, which use rhea as index product and substrate in
    separate columns

    """
    #read file includes rhea id and chEBI id and name
    rhea_dataframe = pd.read_table(file, header=0, index_col=0,
                                          sep='\t')
    #ignore H(+)
    remove_item = ["CHEBI:15378"]
    for id in list(rhea_dataframe.index):
        equation = rhea_dataframe.loc[id,"Equation"]
        subtrates = equation.split("=")[0]
        product = equation.split("=")[1]
        #split by " + ", otherwise H(+) will be sepreate into two part
        subtrate_list = subtrates.split(" + ")
        product_list = product.split(" + ")
        chebiid = (rhea_dataframe.loc[id,"ChEBI identifier"]).split(";")
        participentid = (rhea_dataframe.loc[id,"Participant identifier"]).split(";")
        sub_rh = ""
        sub_ch = ""
        pro_rh = ""
        pro_ch = ""
        for index,item in enumerate(subtrate_list):
            if chebiid[index] not in remove_item:
                sub_rh += participentid[index]+";"
                sub_ch += chebiid[index]+";"
        for i in range(len(product_list),0,-1):
            if chebiid[-i] not in remove_item:
                pro_rh += participentid[-i]+";"
                pro_ch += chebiid[-i]+";"
        rhea_dataframe.loc[id,"substrate_RHEA-COMP"] = sub_rh
        rhea_dataframe.loc[id, "substrate_CHEBI"] = sub_ch
        rhea_dataframe.loc[id,"product_RHEA-COMP"] = pro_rh
        rhea_dataframe.loc[id, "product_CHEBI"] = pro_ch
    return rhea_dataframe

def get_smile(datasets):
    """

    :param datasets: dataframe,which use rhea as index product and substrate in
    separate columns
    :return:
    """

    for i in list(datasets.index):
        substrate = (datasets.loc[i,"substrate_CHEBI"]).split(";")
        #remove the empty string in list
        substrate = list(filter(None, substrate))
        print(substrate)
        products = datasets.loc[i,"product_CHEBI"].split(";")
        products = list(filter(None, products))
        for sub in substrate:
            #some of the mol file is missing in the directory whihc download from rhea
            try:
                mol = Chem.MolFromMolFile("data/mol/"+str(sub.replace(":","_")) +".mol")
            except:
                url = "https://www.ebi.ac.uk/chebi/saveStructure.do?defaultImage=true&chebiId={}&imageId=0".format(sub.split(":")[1])
                # Download from URL.
                with urlopen(url) as webpage:
                    content = webpage.read()
                # Save to file.
                with open("data/mol/"+str(sub.replace(":","_")) +".mol", 'wb') as download:
                    download.write(content)
                print("downloading the missing file: {},pleaze rerun the code".format(sub.replace(":","_")))
            if mol:
                smile = Chem.MolToSmiles(mol)
                mol1 = molecular(sub,smiles=smile)
            print(mol1.get_smiles())

        for pro in products:
            #some of the mol file is missing in the directory whihc download from rhea
            try:
                mol = Chem.MolFromMolFile("data/mol/"+str(pro.replace(":","_")) +".mol")
            except:
                url = "https://www.ebi.ac.uk/chebi/saveStructure.do?defaultImage=true&chebiId={}&imageId=0".format(pro.split(":")[1])
                # Download from URL.
                with urlopen(url) as webpage:
                    content = webpage.read()
                # Save to file.
                with open("data/mol/"+str(pro.replace(":","_")) +".mol", 'wb') as download:
                    download.write(content)
                print("downloading the missing file: {},pleaze rerun the code".format(pro.replace(":","_")))
            if mol:
                smile = Chem.MolToSmiles(mol)
                mol1 = molecular(pro,smiles=smile)

def mol_with_atom_index(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())

    return mol
def main():
    # readfile whihc contains the rhea id and related uniprotid
    #run with command line
    # rh_file = argv[1]
    #run with pycharm
    rh_file = "data/rhea2uniprot_sprot.tsv"
    rheauniprot_dataframe = readrhlist(rh_file)

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

    # df1 = get_substrate_chebi("data/Rhea-ec_2_1_1.tsv")
    # print(df1)
    # get_smile(df1)

    #if the nmbering from rdkit is following the carbon numbering rules?
    mol = mol_with_atom_index('CC1=C[N]C2=C1[C@H](C)C=CC2=O.C(C5=CC4=C3CCN(C3=C(C(=C4[N]5)OC)O)C(=O)N)(=O)N6C7=C(CC6)C8=C(C(=C7O)OC)[N]C(=C8)C(=O)[N]9C%10=C(C=C9)C%11=C(C(=C%10)O)N=CC%11C')
    #Draw.ShowMol(mol,size=(600,600))
    smiles_with_atom_mappings = Chem.MolToSmiles(mol)
    #



    #remove_duplicated_id("data/hmm_out/top_ten_hits_exclude_nomethylrelated/hmmscantbout_top_6_hit.pfam")

    remove_duplicated_id(r"E:\Download\regioselectivity_prediction\data\hmm_out")
    #pikachu.general.draw_smiles(
    #    "CCC(=O)[C@@H]1NC(=O)[C@H](CC(N)=O)NC(=O)C(NC(=O)[C@@H](N)CC(C)C)[C@H](O)c3ccc(Oc2cc1cc(CC)c2O)c(Cl)c3")
if __name__ == "__main__":
    main()
