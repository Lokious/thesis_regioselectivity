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
import subprocess
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
from Bio import AlignIO
from script.sequence import sequences
from script.Model_class import Model_class
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
    remove_item = ["CHEBI:15378"]
    #remove_item = []
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
        j = {}
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
                mol,j = mol_ob.mol_with_atom_index(mol_object=mol,index={})
                smile = Chem.MolToSmiles(mol)
                smiles.append(smile)
                mols_list.append(mol)
            datasets.loc[i,"sub_smiles"]=".".join(smiles)
            datasets.at[i,"sub_mols"] = mols_list
        smiles = []
        mols_list = []
        j = {}
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
                mol,j = mol_ob.mol_with_atom_index(mol_object=mol, index={})
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
    keys = datafrmae_dict.keys()

    return datafrmae_dict

def read_msa_and_encoding(file_name=""):
    """
    simpliest encoding by transfer aminoacide to number and encode it to binary arrary

    :param file:
    :return:
    """
    file = "data/align/Keeplength/{}_align_addmodel".format(file_name)
    import numpy as np
    from sklearn.preprocessing import OneHotEncoder
    align = AlignIO.read(file, "fasta")
    align_array = np.array([list(rec) for rec in align], np.character)
    #the length of aligned sequence
    seq_length = len(align[0])
    print(seq_length)
    char_list = np.unique(align_array)
    char_dictionary = {}
    for i,char in enumerate(char_list):
        char_dictionary[char] = (i+1)


    id = list(rec.id for rec in align)
    align_pd = pd.DataFrame(data=align_array,index=id)
    align_pd = align_pd.drop_duplicates()
    # print(align_pd)

    encode_dictionary = {}
    """"
    for key in char_dictionary.keys():
        align_pd = align_pd.mask(align_pd == key,char_dictionary[key])
        encoding_list = [0]*len(char_dictionary)
        encoding_list[char_dictionary[key]] = 1
        encode_dictionary[char_dictionary[key]] = encoding_list
    """
    # convert number to binary to shorten the size
    for key in char_dictionary.keys():
        align_pd = align_pd.mask(align_pd == key,char_dictionary[key])
        if key != "-":
            bi = [int(x) for x in bin(char_dictionary[key])[2:]]
            encode_dictionary[char_dictionary[key]] = ((6-len(bi))*[0]+bi)
            print(encode_dictionary)
        else:
            encode_dictionary[char_dictionary[key]] = [0]*6

    encode_seq = pd.DataFrame(seq_length * 6 * [],
                              columns=list(range(seq_length * 6)))


    # X = one_hot_encoder.fit_transform(align_pd)
    # print(X)
    #print(len(align_pd.index))
    print(encode_seq.columns)
    for index in align_pd.index:
        line = []
        for aa in list(align_pd.loc[index]):
            line += encode_dictionary[aa]
            print(len(line))
        encode_seq.loc[index] = line

    encode_seq["Entry"] = encode_seq.index
    encode_seq = (encode_seq.reset_index()).drop(columns='index')
    #print(encode_seq)
    with open("data/protein_encoding/{}_simple_encoding_bi".format(file_name),
              "wb") as dill_file:
        dill.dump(encode_seq, dill_file)
    encode_seq.to_csv("data/protein_encoding/{}_simple_encoding_bi.csv".format(file_name))
    return encode_seq

def clean_seq():
    #remove duplicate sequneces in fasta file
    file_list = ["PF05175","PF08241","PF08242","PF13489","PF13649","PF13847"]
    seq = sequences()
    seq.remove_duplicate(file_list)

def merge_uniprot_emebeding():
    file_list = ["PF05175", "PF08241", "PF08242", "PF13489", "PF13649",
                 "PF13847"]
    umiprot_df = pd.DataFrame()
    for file in file_list:
        try:
            with open("data/protein_encoding/{}_simple_encoding_bi".format(file), 'rb') as file1:
                df = dill.load(file1)
            # df = pd.read_csv("data/protein_encoding/{}_simple_encoding.csv".format(file),header=0,index_col=0)
            # print(df)
        except:
            df = read_msa_and_encoding("{}".format(file))
        umiprot_df = pd.concat([umiprot_df,df],axis=0,join='outer')
        #umiprot_df.index=umiprot_df["Entry"]


    "replace the NA with 0, cause the difference in aeq length"
    umiprot_df=(umiprot_df.fillna(int(0)))

    #umiprot_df=umiprot_df.reset_index()
    #some sequences aligned to different hmm,only keep one
    print(umiprot_df.drop_duplicates(subset="Entry",keep="first"))
    umiprot_df=(umiprot_df.drop_duplicates(subset="Entry",keep="first")).reset_index(drop=True)
    umiprot_df.to_csv("data/protein_encoding/protein_seq_simple_encoding_bi.csv")
    return umiprot_df


def create_inputdata(directory, num_bit: int = 2048):
    """

    :param directory:
    :param num_bit:
    :return:
    """
    from script.molecular_class import molecular, reaction
    from script.Model_class import Model_class
    mo_del = Model_class()

    file_name1 = list(glob.iglob(directory + '/seq_smiles_all*'))
    file_name2 = list(glob.iglob(directory + '/diction_atom_all*'))
    print(file_name1)
    print(file_name2)
    assert len(file_name1) == len(file_name2)

    for i, file in enumerate(file_name1):
        name1 = file.split("[")[-1]
        #name2 = (file_name2[i]).split("/")[-1]
        try:
            print("loading data-------")
            X = pd.read_csv("data/input_data/input_dataframe_withoutstructure_dropatoms_{}_{}_withentry_drop_atom_type.csv".format(
                                                      str(num_bit),
                                                      name1), header=0, index_col=0)
            print(X)
        except:
            with open('{}/seq_smiles_all[{}'.format(directory, name1), 'rb') as file1:
                data_with_site = dill.load(file1)
            with open('{}/diction_atom_all[{}'.format(directory, name1), 'rb') as file2:
                diction_atom = dill.load(file2)
            print(diction_atom)
            X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,
                                                      num_bit, 3,
                                                      drop_atoms=True,
                                                      file_name="{}_{}_withentry_drop_atom_type".format(
                                                          str(num_bit),
                                                          name1))
            print(X)
            print(
                "succesfully saved the fingerprint{}".format(name1))
        merge_sequence_and_fingerprint(x=X,
                                       filename=name1,
                                       num_bit=num_bit, add_dataframe="data/protein_encoding/protein_seq_simple_encoding_bi_ri_{}.csv".format(str(num_bit)))
        input_data = pd.read_csv("data/input_data/input{0}fg_dpna_bi_{1}.csv".format(str(num_bit),
                                                            name1), header=0, index_col=0)
        print(input_data)
        # run_model_for_group_data(
        #     input_data, filename="", num_bit=2048)
def merge_sequence_and_fingerprint(x:pd.DataFrame,filename:str="",num_bit:int=2048, add_dataframe=""):

    if add_dataframe == "":
        add_dataframe = pd.read_csv(
            "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv",
            header=0, index_col=0)
    else:
        try:
            add_dataframe = pd.read_csv("{}".format(add_dataframe), header=0, index_col=0)
        except:
            add_dataframe = pd.read_csv(
                "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv",
                header=0, index_col=0)

    print(add_dataframe)
    start_index = num_bit
    #print(start_index)
    if list(add_dataframe.columns)[0]!=str(int(start_index)+1):
        for col in add_dataframe.columns:
            if (col != "Entry") and (col != "index"):
                print(col)
                add_dataframe=add_dataframe.rename(columns={col:str(int(col)+int(start_index)+1)})
            else:
                continue
    add_dataframe.to_csv(
            "data/protein_encoding/protein_seq_simple_encoding_bi_ri_{}.csv".format(str(num_bit)))
    input_dataframe = x.merge(add_dataframe, on="Entry", how="left")

    #drop NA need larger memory
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    col = [i for i in input_dataframe.columns if
           i not in ["Entry", "label", "molecular_id", "methyl_type"]]
    input_dataframe[col] = input_dataframe[col].astype('int32')
    print(input_dataframe)
    input_dataframe = (input_dataframe.reset_index()).drop(columns=["index"])
    print(input_dataframe)
    #
    input_dataframe.to_csv("data/input_data/input{0}fg_dpna_bi_{1}.csv".format(str(num_bit),filename))
    with open("data/input_data/input{0}fg_dpna_bi_{1}".format(str(num_bit),filename), "wb") as dill_file:
        dill.dump(input_dataframe, dill_file)
def run_model_for_group_data(input:pd.DataFrame,filename:str="",num_bit:int=2048):

    mo_del = Model_class()
    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        input)

    mo_del.run_PCA(X_train,y_train,filename)
    X_train = X_train.drop(columns=["methyl_type"])
    X_test = X_test.drop(columns=["methyl_type"])
    y_train = y_train.drop(columns=["methyl_type"])
    y_test = y_test.drop(columns=["methyl_type"])
    model = mo_del.RF_model(X_train, X_test, y_train, y_test,
                            "_{}_{}".format(filename,str(num_bit)))
def main():
    unittest.main()
if __name__ == "__main__":
    main()
