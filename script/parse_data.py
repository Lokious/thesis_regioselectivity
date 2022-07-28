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
import os

import pandas as pd
from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw, AllChem, rdChemReactions
#import pikachu
from script.molecular_class import molecular
import glob
import dill
def target_sequences(file):
    """
    This function is used to return Uniprot entry of target sequences from hmmscan
    :param file:string, file name
    :return: hmmscan_df: dataframe, include all the Uniprot entry from hmmscan output
    """
    name = ['target_name', 'accession_domain','query_name', 'accession',
            'E-value1',  'score1',  'bias1',   'E-value',  'score',  'bias',
            'exp', 'reg', 'clu',  'ov', 'env', 'dom','rep','inc', 'description_of_target']
    hmmscan_df = (pd.read_table(file, header=None, comment='#',names=name,sep= '\s+', skip_blank_lines=True)).dropna()

    hmmscan_df = hmmscan_df.iloc[:,[0,2]]
    for i in hmmscan_df.index:
        hmmscan_df.iloc[i][1] = (hmmscan_df.iloc[i][1]).split("|")[1]
    hmmscan_df.columns = ["index","id"]

    return hmmscan_df
def remove_duplicated_id(directory):
    """
    This is to extract the hits uniprot id from hmmscan and remove the repeat part
    :param directory: sting, directory of hmmscan output
    :return:list of uniprot id
    """
    file_name = glob.iglob(directory+'/*.tsv')
    datafrmae_list = []
    for i in file_name:
        df = target_sequences(i)
        datafrmae_list += df['id'].tolist()
        #return uniprot id
        return list(set(datafrmae_list))

def readrhlist(file_name):

    rheauniprot_dataframe = pd.read_table(file_name, header=0, index_col=3,sep='\t', encoding="ISO-8859-1")
    rheauniprot_dataframe = rheauniprot_dataframe.rename_axis(index='Entry')
    rheauniprot_dataframe["RHEA_ID"] = rheauniprot_dataframe["RHEA_ID"].apply(lambda x: str(x))
    return rheauniprot_dataframe


def read_sequence(file):
    lines = pd.read_excel(file,index_col=0)
    return lines


def rheaid_touniprot(rhea,data_frame):

    entrys = pd.DataFrame()
    for index,id in enumerate(rhea):
        print(id)
        number = (id.strip()).split(":")[-1]
        uniprot_entry = data_frame[data_frame['RHEA_ID'] == int(number)]
        print(uniprot_entry)
        entrys.loc[index, 'RHEA_ID'] = id
        entrys.loc[index, "Entry"] = uniprot_entry
        if uniprot_entry.empty:
            entrys.loc[index,'RHEA_ID'] = id
            entrys.loc[index,"Entry"] = 'NA'
        else:
            continue
    print(entrys)
    return entrys
def rheaid_to_uniprot(uniprotids,data_frame):

    entrys = pd.DataFrame()
    entrys["Entry"] = uniprotids
    entrys = data_frame.join(entrys)
    entrys["Entry"] = entrys.index
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
    #remove_item = ["CHEBI:15378"]
    remove_item = []
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

def get_smile(datasets)->pd.DataFrame:
    """

    :param datasets: dataframe,which use rhea as index product and substrate in
    separate columns
    :return: datdframe which include smiles and mol objects
    """
    datasets["sub_mols"] = pd.DataFrame(
        len(datasets["substrate_CHEBI"]) * [0]).astype('object')
    datasets["pro_mols"] = pd.DataFrame(
        len(datasets["substrate_CHEBI"]) * [0]).astype('object')
    for i in list(datasets.index):
        substrate = (datasets.loc[i,"substrate_CHEBI"]).split(";")
        #remove the empty string in list
        substrate = list(filter(None, substrate))
        products = datasets.loc[i,"product_CHEBI"].split(";")
        products = list(filter(None, products))
        smiles = []
        mols_list = []
        j = 0
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
                mol_ob = molecular()
                mol,j = mol_ob.mol_with_atom_index(mol_object=mol,index=j)
                smile = Chem.MolToSmiles(mol)
                smiles.append(smile)
                mols_list.append(mol)
            datasets.loc[i,"sub_smiles"]=".".join(smiles)
            datasets.at[i,"sub_mols"] = mols_list
        smiles = []
        mols_list = []
        j = 0
        for pro in products:
            #some of the mol file is missing in the directory which download from rhea

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
                mol_ob = molecular()
                mol,j = mol_ob.mol_with_atom_index(mol_object=mol, index=j)
                smile = Chem.MolToSmiles(mol)
                smiles.append(smile)
                mols_list.append(mol)
            datasets.loc[i,"pro_smiles"] = ".".join(smiles)
            datasets.at[i,"pro_mols"] = mols_list

    return datasets


def merge_reaction(dataset)->pd.DataFrame:
    """
    add the rnx to dataframe
    :return: pd.DataFrame
    """
    dataset["rxn"] = pd.DataFrame(
        len(dataset["RHEA_ID"]) * [0]).astype('object')
    dataset["rxn_smart"] = pd.DataFrame(
        len(dataset["RHEA_ID"]) * [0]).astype('object')
    #print(dataset)

    for i, rheaid in enumerate(dataset["RHEA_ID"]):
        # read rxn from smile
        reaction_smile ="{}>>{}".format(dataset.loc[i,"sub_smiles"],dataset.loc[i,"pro_smiles"])
        rxn = AllChem.ReactionFromSmarts(reaction_smile, useSmiles=True)
        if rxn != None:
            dataset.loc[i,"rxn"] = rxn
            dataset.loc[i,"rxn_smart"] = rdChemReactions.ReactionToSmarts(rxn)
        else:
            dataset.loc[i, "rxn"] = None
    return dataset

def read_mannual_data(file=r"E:\Download\regioselectivity_prediction\data\mannual_data.csv"):
    """
    read manually made data and save the structure with atom index
    :param file:
    :return:
    """
    manual_data=pd.read_csv(file,header=0,index_col=0,encoding="ISO-8859-1")
    manual_data["index_smile"] = pd.DataFrame(
        len(manual_data.index) * [0]).astype('object')
    manual_data["mainsub_mol"] = pd.DataFrame(
        len(manual_data.index) * [0]).astype('object')
    methyl_site_dictionary = {}
    for index in manual_data.index:
        smile = manual_data.loc[index,'substrate_smiles']
        mol_ob = molecular()
        mol, j = mol_ob.mol_with_atom_index(smile=smile)
        try:
            methyl_site_dictionary[index] = (manual_data.loc[index,'reactant_site']).split(",")
        except:
            methyl_site_dictionary[index]=[]
        if mol:
            smile = Chem.MolToSmiles(mol)
            manual_data.loc[index, "mainsub_mol"] = mol
            manual_data.loc[index, "index_smile"] = smile
            img=Draw.MolToImage(mol,size=(600,600))
            file_name = "data/manual/{}_{}.png".format(manual_data.loc[index,'Entry'],index)
            img.save(file_name)
        else:
            manual_data.loc[index, "mainsub_mol"] = "NA"
            manual_data.loc[index, "index_smile"] = smile
    #save panda.Dataframe object
    with open("data/mannual_data", "wb") as dill_file:
        dill.dump(manual_data, dill_file)
    #save csv file to for check
    manual_data.to_csv("data/mannual_data.csv")
    with open("data/methyl_site_dictionary", "wb") as dill_file:
        dill.dump(methyl_site_dictionary, dill_file)

def group_by_domain(directory,dataframe):

    file_name = glob.iglob(directory + '/*.tsv')
    datafrmae_dict = {}
    for i in file_name:

        df = target_sequences(i)
        seed = (i.split("\\")[-1]).split("_")[0]
        datafrmae_dict[seed] = dataframe.loc[dataframe["Entry"].isin(list(df["id"]))]
    with open("data/dictionary_with_sepreate_datasets_by_seeds", "wb") as dill_file:
        dill.dump(datafrmae_dict, dill_file)

    for key in datafrmae_dict.keys():
        f = open("data/{}.fasta".format(key), "a")
        for index in datafrmae_dict[key].index:

            entry = datafrmae_dict[key].loc[index,"Entry"]
            f.write(">{}\n".format(entry))
            f.write("{}\n".format(datafrmae_dict[key].loc[index,"Sequence"]))
        f.close()
        run_maff(datafrmae_dict.keys())
    return datafrmae_dict

def run_maff(files):
    for file in files:
        cmd1 = "nohup mafft - -auto - -anysymbol - -add {0}.fasta - -reorder - -thread - 2 hmm_out/top_ten_hits_exclude_nomethylrelated/{0}_seed.fasta_aln.txt > align/{0}_align_addmodel &\n".format(file)
        os.system(cmd1)
        cmd2=" nohup mafft --auto --anysymbol --addfragments {0}.fasta --reorder --thread -2 hmm_out/top_ten_hits_exclude_nomethylrelated/{0}_seed.fasta_aln.txt > align/{0}_align &\n".format(file)
        os.system(cmd2)
def main():
    unittest.main()
if __name__ == "__main__":
    main()
