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
import typing
from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw, AllChem, rdChemReactions
#import pikachu
from upsetplot import plot
from molecular_class import Molecule
import glob
import dill
from Bio import AlignIO, SeqIO
from sequence import Sequences
import numpy as np
import unittest
import  sys

def target_sequences(file):
    """
    This function is used to return Uniprot entry of target sequences from hmmscan

    :param file:string, file name
    :return: hmmscan_df: dataframe, include all the Uniprot entry from hmmscan output
    """
    name = ['target_name', 'accession_domain','query_name', 'accession',
            'E-value1',  'score1',  'bias1',   'E-value',  'score',  'bias',
            'exp', 'reg', 'clu',  'ov', 'env', 'dom','rep','inc', 'description_of_target']
    #hmmscan_df = (pd.read_table(file, header=None, comment='#',names=name,sep= '\s+', skip_blank_lines=True)).dropna()
    hmmscan_df=pd.read_csv(file,sep="\t")
    hmmscan_df = hmmscan_df.iloc[:,[0,2]]
    for i in hmmscan_df.index:
        hmmscan_df.iloc[i][1] = (hmmscan_df.iloc[i][1]).split("|")[1]
    hmmscan_df.columns = ["index","id"]

    return hmmscan_df

def get_fasta_file_from_hmmsearch_hit( input: typing.Optional[typing.Union[str, pd.DataFrame]] = None,domain=""):
    if isinstance(input, str):
        hmmsearch_df = read_hmmsearch_out(input)
        print(hmmsearch_df)
        hmmsearch_df.drop_duplicates(subset=['entry'], inplace=True)
        seq_entry_df = pd.read_csv("../autodata/rawdata/uniprot-ec2.1.1.tsv",header=0,sep="\t")
        print(seq_entry_df)
        seq_entry_df = seq_entry_df.loc[:,["Entry","Sequence"]]
        seq_entry_df.index = seq_entry_df["Entry"]
        seq_entry_df.drop(columns="Entry",inplace=True)
        print(seq_entry_df)
        entry_list=[]
        file = open("../autodata/sequences/{}.fasta".format(domain), "w")
        for i in hmmsearch_df.index:
            entry = hmmsearch_df.loc[i,"entry"]
            if entry not in entry_list:
                file.write(">{}\n".format(entry))
                file.write("{}\n".format(seq_entry_df.loc[entry, "Sequence"]))
        else:
            print("File saved.\n ../autodata/sequences/{}.fasta".format(domain))
            file.close()
    elif isinstance(input, pd.DataFrame):
        print("input is dataframe")
        print(input)
        seq_entry_df = pd.read_csv("../autodata/rawdata/uniprot-ec2.1.1.tsv",
                                   header=0, sep="\t")
        seq_entry_df = seq_entry_df.loc[:, ["Entry", "Sequence"]]
        seq_entry_df.index = seq_entry_df["Entry"]
        seq_entry_df.drop(columns="Entry", inplace=True)
        print(seq_entry_df)
        entry_list = []
        file = open("../autodata/sequences/{}.fasta".format(domain), "w")
        for i in input.index:
            entry = input.loc[i, "entry"]
            if entry not in entry_list:
                file.write(">{}\n".format(entry))
                file.write("{}\n".format(seq_entry_df.loc[entry, "Sequence"]))
        else:
            file.close()
def read_hmmsearch_out(file):
    """

    :param file: hmmsearch domout.tsv
    :return:hmmsearch_df
    """
    print(file)
    file = open(file,encoding='utf-8').readlines()

    hmmsearch_df={}
    hmmsearch_df["domain"]=[]
    hmmsearch_df["entry"] = []
    print(len(hmmsearch_df.keys()))
    for i,line in enumerate(file):
        if line.startswith("#")==False:
            #print(line.split())
            entry = line.split()[0]
            domain = line.split()[4]
            hmmsearch_df["entry"].append(entry.split("|")[1])
            hmmsearch_df["domain"].append(domain)
    hmmsearch_df=pd.DataFrame(hmmsearch_df,index=range(len(hmmsearch_df["domain"])))
    #print((hmmsearch_df["domain"].value_counts()[:10]))
    return hmmsearch_df

def upsetplot(seq_domain_df,version):
    """
    This is the function for plot the upsetplot for domains and sequences overlap
    :param seq_domain_df: pd.Dataframe includes two columns["entry","domain"]
    :return: None
    """
    #only  show the top 15 domains
    columns=list((seq_domain_df["domain"].value_counts()[:15]).index)
    print(columns)
    new_df = pd.DataFrame(columns=columns, index=seq_domain_df["entry"])
    print(new_df)
    for index_old in seq_domain_df.index:
        entry=seq_domain_df.loc[index_old,"entry"]
        domain = seq_domain_df.loc[index_old, "domain"]
        if domain in columns:
            new_df.loc[entry,domain] = True
    else:
        new_df.fillna(False,inplace=True)
        print(new_df)
        #new_df.to_csv("new_df.csv")
    from upsetplot import plot as upset_plot
    from upsetplot import from_indicators
    #new_df=pd.read_csv("new_df.csv",header=0,index_col=0)
    print(new_df)
    input_plotdf = from_indicators(new_df)
    input_plotdf.sort_values(ascending=False, inplace=True)
    print(input_plotdf)

    upset_plot(input_plotdf, subset_size='count',sort_by='cardinality',show_counts=True)
    from matplotlib import pyplot
    current_figure = pyplot.gcf()
    current_figure.savefig("upset_{}.png".format(version))

def extract_pdb_structure_from_hmmer_output(domain="PF13847",path="../autodata/align/separate_by_domain/"):
    """
    This function is use to find the pdb structure include active site and related to methyl
    :param domain:
    :param path:
    :return:
    """
    file=open("{}pdb_{}_trim.hmmer".format(path,domain)).readlines()

    pdb_df={}
    pdb_df["pdb_entry"]=[]
    pdb_df["description"] = []
    print(file)
    start=False
    for i,line in enumerate(file):
        if line.startswith("#")==False:
            if line == "\n":
                print("emptyline")
            if start==True:
                if line == "\n":
                    pdb_df = pd.DataFrame(pdb_df, index=range(
                        len(pdb_df["pdb_entry"])))
                    print(pdb_df)
                    return pdb_df
                #print(line.split())
                pdb_entry = line.split()[8]
                description = "".join(line.split()[9:])
                pdb_df["pdb_entry"].append(pdb_entry)
                pdb_df["description"].append(description)
            if line.startswith("-")==True:
                start=True
    pdb_df=pd.DataFrame(pdb_df,index=range(len(pdb_df["pdb_entry"])))
    print(pdb_df)

def sepreate_sequence_based_on_domain_without_overlap(seq_domain_df):
    """
    This function is use to extract sequences only has one domian for top 10 domains

    :param seq_domain_df: pd.Dataframe with domain and uniprot entry
    :return: None
    """
    print(seq_domain_df["domain"].value_counts()[:15])
    seq_domain_df.drop_duplicates(subset=['entry'], inplace=True)
    print(seq_domain_df["domain"].value_counts()[:15])
    print(seq_domain_df)
    domains = list(seq_domain_df["domain"].value_counts()[:10].index)
    print(domains)

    seq_entry_df=pd.read_csv("../autodata/rawdata/uniprot-ec2.1.1.tsv",header=0,sep="\t")
    print(seq_entry_df)
    seq_entry_df=seq_entry_df.loc[:,["Entry","Sequence"]]
    seq_entry_df.index=seq_entry_df["Entry"]
    seq_entry_df.drop(columns="Entry",inplace=True)
    print(seq_entry_df)
    for domain in domains:
        file = open("../autodata/sequences/no_overlap_sequences/{}.fasta".format(domain), "w")
        hmm_df=seq_domain_df.loc[seq_domain_df['domain'] == domain]
        print(hmm_df)
        for i in hmm_df.index:
            entry = hmm_df.loc[i,"entry"]
            file.write(">{}\n".format(entry))
            file.write("{}\n".format(seq_entry_df.loc[entry, "Sequence"]))
        else:
            file.close()
    return domains

def save_sequences_from_hmmscan_result(file="../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/all_sequences_build_hmm_search.tsv"):
    file = open(file, encoding='utf-8').readlines()

    hmmscan_df = {}

    hmmscan_df["entry"] = []
    for i, line in enumerate(file):
        if line.startswith("#") == False:

            entry = line.split()[3]

            hmmscan_df["entry"].append(entry.split("|")[1])
    hmmscan_df = pd.DataFrame(hmmscan_df,
                                index=range(len(hmmscan_df["entry"])))
    # print((hmmscan_df["domain"].value_counts()[:30]))
    get_fasta_file_from_hmmsearch_hit(hmmscan_df,"all_sequences_build_hmm_search")
    hmmscan_df.drop_duplicates(subset=["entry"],inplace=True)
    print(hmmscan_df)
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
                mol = Chem.MolFromMolFile("../data/mol/"+str(sub.replace(":","_")) +".mol")
            except:
                url = "https://www.ebi.ac.uk/chebi/saveStructure.do?defaultImage=true&chebiId={}&imageId=0".format(sub.split(":")[1])
                # Download from URL.
                with urlopen(url) as webpage:
                    content = webpage.read()
                # Save to file.
                with open("../data/mol/"+str(sub.replace(":","_")) +".mol", 'wb') as download:
                    download.write(content)
                print("downloading the missing file: {},pleaze rerun the code".format(sub.replace(":","_")))
            if mol:
                mol_ob = Molecule()
                mol,j = mol_ob.mol_with_atom_index(mol,index={})
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
                mol = Chem.MolFromMolFile("../data/mol/"+str(pro.replace(":","_")) +".mol")
            except:
                url = "https://www.ebi.ac.uk/chebi/saveStructure.do?defaultImage=true&chebiId={}&imageId=0".format(pro.split(":")[1])
                # Download from URL.
                with urlopen(url) as webpage:
                    content = webpage.read()
                # Save to file.
                with open("../data/mol/"+str(pro.replace(":","_")) +".mol", 'wb') as download:
                    download.write(content)
                print("downloading the missing file: {},pleaze rerun the code".format(pro.replace(":","_")))
            if mol:
                mol_ob = Molecule()
                mol,j = mol_ob.mol_with_atom_index(mol, index={})
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
        mol_ob = Molecule()
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


def read_msa_and_encoding(file_name=""):
    """
    simpliest encoding by transfer aminoacide to number and encode it to binary arrary

    :param file_name: string, the file name for aligned sequences
    :return:
    """
    file = "../autodata/align/Keeplength/{}_addmodel_rm".format(file_name)

    #read alignment
    #align = AlignIO.read(file_name, "clustal")
    align = AlignIO.read(file, format="fasta")
    align_array = np.array([list(rec) for rec in align], np.character)
    print(align_array)
    id = list(rec.id for rec in align)
    align_pd = pd.DataFrame(data=align_array,index=id)
    align_pd = align_pd.drop_duplicates()

    char_list = np.unique(align_array)

    char_dictionary = {}
    for i,char in enumerate(char_list):
        char_dictionary[char] = (i+1)




    #drop columns which for over 80% sequences is gap
    from collections import Counter
    for column in list(align_pd.columns):
        gap_num = Counter(align_pd[column])[b"-"]
        gap_percentage = (gap_num/len(align_pd.index))

        if gap_percentage > 0.8:
            align_pd.drop(columns=[column],inplace=True)
    print(align_pd)
    # the length of aligned sequence after remove some columns contains a lot of gaps
    seq_length=len(align_pd.columns)
    print("seq_length:{}".format(seq_length))

    encode_dictionary = {}

    for key in char_dictionary.keys():
        align_pd = align_pd.mask(align_pd == key,char_dictionary[key])
        encoding_list = [0]*len(char_dictionary.keys())
        encoding_list[char_dictionary[key]-1] = 1
        encode_dictionary[char_dictionary[key]] = encoding_list
    #print(encode_dictionary)
    encode_seq = pd.DataFrame(seq_length * len(char_list) * [],
                              columns=list(range(seq_length * len(char_list))))


    for index in align_pd.index:
        line = []
        for aa in list(align_pd.loc[index]):
            line += encode_dictionary[aa]
        encode_seq.loc[index] = line

    encode_seq["Entry"] = encode_seq.index
    # encode_seq["domain"] = len(encode_seq.index)*[file_name]
    encode_seq = (encode_seq.reset_index()).drop(columns='index')
    print(encode_seq)
    # with open("../autodata/protein_encoding/{}_onehot_encoding_sepreate".format(file_name),
    #           "wb") as dill_file:
    #     dill.dump(encode_seq, dill_file)
    encode_seq.to_csv("../autodata/protein_encoding/{}_onehot_encoding_sepreate.csv".format(file_name))
    # encode_seq.to_csv(
    #     "{}_onehot_encoding_sepreate_without_drop.csv".format(
    #         file_name))
    return encode_seq
def K_mer_generation(k:int=0)->list:
    """

    :param k:
    :return:
    """
    from itertools import product
    AA=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S","U", "T", "W", "Y", "V","B","Z","X","J","O","-"]
    k_mer = list(product(AA,repeat=k))
    for i,item in enumerate(k_mer):
       k_mer[i]="".join(item)
    return k_mer
def k_mer_encoding(file_name,k):
    """
    This is the function of create k_mer_encoding for sepreate groups
    :param file_name: string,group name
    :param k:int, the length of k
    :return: align_pd:sequences encoding pd.Dataframe with group and id
    """
    file = "../autodata/sequences/{}_rm.fasta".format(file_name)
    # read alignment
    seqs = SeqIO.parse(file, "fasta")
    sequences_dictionary={}
    for record in seqs:
        sequences_dictionary[record.id]=("".join(list(record.seq)))
    #print(sequences_dictionary)
    #all possiable k_mer
    K_mer_list=K_mer_generation(k)
    ids = list(sequences_dictionary.keys())
    #print("id:{}".format(ids))
    align_pd = pd.DataFrame(data=np.zeros((len(ids),len(K_mer_list+["ID","group"]))),columns=(K_mer_list+["ID","group"]),index=list(set(ids)))
    align_pd["group"]=align_pd["group"].astype("string")
    #print(align_pd)
    for index,id in enumerate(list(align_pd.index)):
        sequences = sequences_dictionary[id]
        print(index)
        for start_site in range(0,len(sequences),k):
            three_mer= "".join(sequences[start_site:(start_site + k)])
            if len(three_mer)== k:
                align_pd.loc[id,three_mer]=int(align_pd.loc[id,three_mer]+1)
                align_pd.loc[id,"group"]=file_name
            else:
                continue

    align_pd["Entry"] = align_pd.index
    align_pd = (align_pd.reset_index()).drop(columns='index')
    print(align_pd)
    print(align_pd.dtypes)
    align_pd.to_csv("../autodata/protein_encoding/{}_k_mer_encoding_without_align_26_08.csv".format(file_name))
    return align_pd

def use_atom_properties_for_sequences_encoding(input_df:pd.DataFrame=None,file_name:str=None,structure_chain="3rod.pdb_chainA_s001",start:int=0,group:str="",file_format="fasta",pdb_name="3rod.pdb"):
    #add properties as features
    seq = Sequences()
    if input_df!=None:
        print(input_df)
    else:
        try:
            print("use input dataframe")
            input_df=pd.read_csv("../autodata/protein_encoding/active_site/{}_active_site_df.csv".format(group),header=0,index_col=0)
            print(input_df)
        except:
            input_df = seq.get_sites_from_alignment(fileformat=file_format,pdb_name=pdb_name,
                file=file_name,structure_chain=structure_chain,
                start_pos=start,group=group)
            print(input_df)


    properties=['charge', 'volume', 'hydropathy', 'hydrophobicity','similarity_score']
    # create empty dataframe
    property_df=pd.DataFrame(index=input_df.index)
    for col in input_df.columns:
        for pro in properties:
            col_name = str(col)+pro
            property_df[col_name] = ["NA"]*len(property_df.index)

    # calculate properties and save to dataframe
    for i in property_df.index:

        for col in input_df.columns:
            structure_aa = input_df.loc[structure_chain, col]
            aa = input_df.loc[i, col]
            propertites_dictionary=seq.amino_acid_properties(amino_acid=aa, structure_aa=structure_aa)
            for pro in properties:
                col_name = str(col) + pro
                property_df.loc[i, col_name] = propertites_dictionary[pro]
    print(property_df)
    property_df.to_csv("../autodata/protein_encoding/active_site/{}_AA_properties_encoding.csv".format(group))
    return property_df

def merge_encoding_to_onefile():

    files = ["O","N","S","C","Co"]
    pd_list = []
    for file in files:
        try:
            p = pd.read_table("jdb")
            #pd_add = pd.read_csv("../autodata/protein_encoding/{}_k_mer_encoding_without_align.csv".format(file),header=0,index_col=0)
            #print(pd_add)
        except:
            print("create k_mer encoding.....")
            pd_add=k_mer_encoding(file, 3)

        pd_list.append(pd_add.astype(object))
    all_pd=pd.concat(pd_list,axis=0)
    print(all_pd)
    print(len(all_pd.sum() != 0.0))
    print(all_pd.columns[0])
    all_pd=(all_pd.loc[:, (all_pd.sum() != 0.0)])
    #print(all_pd)
    all_pd.to_csv("../autodata/protein_encoding/all_k_mer_encoding_sepreate_without_align.csv")
    return all_pd

def get_active_and_binding_site_realted_to_methylation():
    """
    The input data is downloaded from Uniprot database
    :return:
    """
    entry_acs_pd = pd.read_excel("../autodata/uniprot_activesite.xlsx",header=0,index_col=None)
    entry_acs_pd["active_site_AA"]=pd.DataFrame(
            len(entry_acs_pd.index) * [0])
    entry_acs_pd["binding_site_AA"]=pd.DataFrame(
            len(entry_acs_pd.index) * [0])
    entry_acs_pd["CHEBI"]=pd.DataFrame(
            len(entry_acs_pd.index) * [0])
    entry_acs_pd["PDB_binding_site"]=pd.DataFrame(
            len(entry_acs_pd.index) * [0])
    for index in entry_acs_pd.index:
        active_sites_list=(entry_acs_pd.loc[index,"Active site"]).split("ACT_SITE")[1:]
        siteAA_list=[]
        #split active site and extract related to methylation
        for active_site_list in active_sites_list:
            site = active_site_list.split(";")[0]

            if ("methyl" or "MT") in "".join(active_site_list.split(";")):

                siteAA_list.append(site)
        else:
            entry_acs_pd.loc[index,"active_site_AA"]=",".join(siteAA_list)
        # split binding site and save ligand's CHEBI, save intwo columns in same order

        if pd.isna(entry_acs_pd.loc[index,"Binding site"]) == False:

            binding_sites_list = (entry_acs_pd.loc[index,"Binding site"]).split("BINDING")[1:]
            bindingsiteAA_list = []
            chebi_list = []
            pdb_id_list = []
            print(index)

            for binding_site_list in binding_sites_list:

                site = binding_site_list.split(";")[0]
                bindingsiteAA_list.append(site)
                try:
                    chebi_list.append((binding_site_list.split('"ChEBI:')[1]).split(";")[0])
                except:
                    chebi_list.append(" ")
                try:
                    pdb_id_list.append((binding_site_list.split("PDB:")[1])[:4])
                    #print(binding_site_list.split("PDB:")[1])
                    print(pdb_id_list)
                except:
                    pdb_id_list.append(" ")
            entry_acs_pd.loc[index, "binding_site_AA"]=",".join(bindingsiteAA_list)
            entry_acs_pd.loc[index, "CHEBI"] = ",".join(chebi_list)
            entry_acs_pd.loc[index, "PDB_binding_site"] = ",".join(pdb_id_list)
        else:
            entry_acs_pd.loc[index, "binding_site_AA"] = ""
            entry_acs_pd.loc[index, "CHEBI"] = "NA"
            entry_acs_pd.loc[index, "PDB_binding_site"] = ""
    #drop those without related active sites
    entry_acs_pd=entry_acs_pd[~((entry_acs_pd["active_site_AA"] == "") & (entry_acs_pd["binding_site_AA"] == ""))]
    #entry_acs_pd.dropna(subset=["active_site_AA","binding_site_AA"],inplace=True)
    print(entry_acs_pd)
    entry_acs_pd.to_csv("../autodata/entry_with_activesite.csv")

def merge_active_site_and_methyltype (activesite_file, fingerprint_file):
    active_site_df= pd.read_csv(activesite_file, header=0, index_col=0)
    #active_site_df=active_site_df.fillna(0)
    active_site_df=pd.DataFrame(active_site_df,columns=["active_site_AA", "Entry"])
    fg_df=pd.read_csv(fingerprint_file,header=0,index_col=0)
    merg_df=fg_df.merge(active_site_df, on="Entry", how="left")
    print(merg_df)
    merg_df.dropna(inplace=True)
    print(merg_df)
    merg_df.to_csv("../autodata/AC_merge.csv")

def clean_seq():
    #remove duplicate sequneces in fasta file
    file_list = ["PF08241", "PF05175", "PF08242", "PF13489", "PF13649",
                 "PF13847"]
    seq = Sequences()
    seq.remove_duplicate(file_list)

def read_fasta_file(file_name=""):
    """
    This function is reading fasta file into dictionary

    :param file_name: path and name of a fasta file
    :return: dictionary of sequences
    """
    file = open(file_name).readlines()
    seq_dictionary = {}
    seq_name = ""
    for line in file:
        if line.startswith(">"):
            seq_name = "".join(list(line.strip("\n"))[1:])
            seq_dictionary[seq_name]=""
            continue
        seq_dictionary[seq_name] += line.strip("\n")
    return seq_dictionary

def check_sequences_similarity(fasta_file=""):
    """
    This function is to compare the sequences similarity
    :return:
    """
    import difflib
    groups=["O","S","N","C"]
    for group in groups:
        seq_dictionary=read_fasta_file("../autodata/align/Keeplength/{}_seed_addmodel_rm".format(group))
        similarity_matrix=pd.DataFrame(index=list(seq_dictionary.keys()),columns=list(seq_dictionary.keys()))
        visit_list=[]
        for column in similarity_matrix.columns:
            for index in similarity_matrix.index:
                if {column,index} not in visit_list:
                    seq1=seq_dictionary[column]
                    seq2=seq_dictionary[index]
                    visit_list.append({column,index})
                    similarity_matrix.loc[index,column]=difflib.SequenceMatcher(None,seq1,seq2 ).ratio()
                    print(similarity_matrix)
                else:
                    similarity_matrix.loc[index, column]=similarity_matrix.loc[column,index]
        print(similarity_matrix)
        similarity_matrix.to_csv("../autodata/sequences/{}_rm_similarity_matrix.csv".format(group))

# def merge_uniprot_emebeding():
#     file_list = ["PF08241","PF05175",  "PF08242", "PF13489", "PF13649",
#                  "PF13847"]
#     umiprot_df = pd.DataFrame()
#     for file in file_list:
#         try:
#             with open("data/protein_encoding/{}_onehot_encoding".format(file), 'rb') as file1:
#                 df = dill.load(file1)
#             # df = pd.read_csv("data/protein_encoding/{}_simple_encoding.csv".format(file),header=0,index_col=0)
#             # print(df)
#         except:
#             df = read_msa_and_encoding("{}".format(file))
#         umiprot_df = pd.concat([umiprot_df,df],axis=0,join='outer')
#         #umiprot_df.index=umiprot_df["Entry"]
#
#
#     "replace the NA with 0, cause the difference in aeq length"
#     umiprot_df=(umiprot_df.fillna(int(0)))
#
#     #umiprot_df=umiprot_df.reset_index()
#     #some sequences aligned to different hmm,only keep one
#     print(umiprot_df.drop_duplicates(subset="Entry",keep="first"))
#     umiprot_df=(umiprot_df.drop_duplicates(subset="Entry",keep="first")).reset_index(drop=True)
#     umiprot_df.to_csv("data/protein_encoding/protein_seq_simple_encoding_bi.csv")
#     return umiprot_df



def main():
    use_atom_properties_for_sequences_encoding(file_name="../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/PF08241.15PF03602.18/PF08241.15PF03602.18_hmmalign_out_pdb_5WP4.aln",group="PF08241.15PF03602.18",file_format="clustal",start=0, structure_chain="5WP4_1|Chain",pdb_name="5wp4.pdb")
    #unittest.main()
    #, min_subset_size = 100, max_degree = 4

    #upsetplot(seq_domain_df,"ec211_represent_seq")
    # versions = ["Pfam25.0", "Pfam29.0", "Pfam35.0"]
    # for version in versions:
    #     seq_domain_df=read_hmmsearch_out("../autodata/align/{}uniprot_2_1_1_domout.tsv".format(version))
    #     upsetplot(seq_domain_df,version)

    # extract_pdb_structure_from_hmmer_output(domain="PF13847",
    #                                         path="../autodata/align/separate_by_domain/")

    #check_sequences_similarity()

    # use_atom_properties_for_sequences_encoding(file_name="../autodata/align/align_seed_sequences_with_structure/C_4u1q_align_sequences",structure_chain="4u1q.pdb_chainA_s001",
    #     start=4,group="C")
    #merge_active_site_and_methyltype("../autodata/entry_with_activesite.csv","../autodata/fingerprint_bit128_radius3_all_data_drop_atom.csv")
    # read_msa_and_encoding(file_name="N_seed")
    #merge_encoding_to_onefile()
    # for i in ["O","N","O_N","S","C","Se","Co"]:
    #      trget_df=read_hmmscan_out("../autodata/align/hmmsearch/{}_tblout.tsv".format(i))
    # get_active_and_binding_site_realted_to_methylation()
    #merge_encoding_to_onefile()
if __name__ == "__main__":
    main()
