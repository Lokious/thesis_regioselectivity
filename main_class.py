#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
datafile should be put in the directory data/ ,this script should be under regioselectivity_prediction
includes Rhea-ec_2_1_1.tsv
This is the main function
"""
import dill

from script import parse_data
from script.Model_class import Model_class
import pandas as pd
import numpy as np
import copy
import time
from datetime import date
from rdkit.Chem import Draw,DataStructs
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from rdkit import Chem
import glob
from script.sequence import sequences


def run_model_for_group_data(input:pd.DataFrame,filename:str="",num_bit:int=2048):

    mo_del = Model_class()
    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        input)

    #mo_del.run_PCA(X_train,y_train,filename)
    X_train = X_train.drop(columns=["methyl_type"])
    X_test = X_test.drop(columns=["methyl_type"])
    y_train = y_train.drop(columns=["methyl_type"])
    y_test = y_test.drop(columns=["methyl_type"])
    model1 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                            "_{}_{}".format(filename,str(num_bit)))
    # model2 = mo_del.SVM(X_train, X_test, y_train, y_test,
    #                         "_{}_{}".format(filename,str(num_bit)))

def build_different_input(x="",num_bit:int=0,radius:int=0,seqfile:str="6_seed_onehot_encoding.csv",group=""):
    """

    :param num_bit:
    :param radius:
    :return:
    """
    mo_del = Model_class()

    try:
        if x=="":
            X = pd.read_csv(
                "data/input_dataframe_withoutstructure_dropatoms{}_drop_duplicate_drop_atom_withtype_bond{}.csv".format(
                    str(num_bit), str(radius)), header=0, index_col=0)
        else:
            X = pd.read_csv("{}_{}_withtype_bond{}_['{}'].csv".format(x,num_bit,radius,group),header=0,index_col=0)
    except:
        data_with_site = pd.read_csv("data/seq_smiles_all_MANUAL.csv", header=0, index_col=0)
        with open('data/diction_atom_all', 'rb') as file1:
            diction_atom = dill.load(file1)
        X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,
                                                  num_bit, radius,
                                                  drop_atoms=True,
                                                  file_name="{}_drop_duplicate_drop_atom_withtype_bond{}".format(str(num_bit),str(radius)))
    print(X)

    add_dataframe = pd.read_csv("data/protein_encoding/{}".format(seqfile), header=0, index_col=0)
    print(add_dataframe)

    start_index = num_bit*2
    #print(start_index)
    if list(add_dataframe.columns)[0] != str(str):
        print("renaming columns----")
        map_dictionary ={}
        for col in add_dataframe.columns:
            if (col != "Entry") and (col != "index"):

                map_dictionary[col] = str(int(col)+int(start_index))

            else:
                continue
        add_dataframe = add_dataframe.rename(columns=map_dictionary)
        print("####rename column finished####")
        print(add_dataframe)
        add_dataframe.to_csv(
                "data/protein_encoding/{}_{}fg.csv".format(seqfile,str(num_bit)))
    print("merging fingerprint and sequences encoding---------")
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    print(input_dataframe)
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    print(input_dataframe)
    print("saving input data.......")
    #
    input_dataframe.to_csv("data/input_data/input{}fg_dpna_bond{}_{}.csv".format(str(num_bit),str(radius),seqfile))

def sepreate_input(file="",numbit:int=2048,bond:int=2):
    input_dataall=pd.read_csv(file,header=0,index_col=0)
    seprate_dataset = input_dataall.groupby(by=["methyl_type"])
    for group in seprate_dataset.groups:
        sub_df = seprate_dataset.get_group(group)
        group = sub_df["methyl_type"].unique()
        print(group)
        sub_df.reset_index(drop=True,inplace=True)

        sub_df.to_csv("data/input_data/group/{}_{}_{}.csv".format(group,str(numbit),str(bond)))
def main():
    today = date.today()
    # dd/mm/YY
    d1 = today.strftime("%d_%m_%Y")
    mo_del = Model_class()
    data_with_site = pd.read_csv("data/mannual_data.csv", header=0,
                                 index_col=0)
    with open('data/methyl_site_dictionary', 'rb') as file1:
        diction_atom = dill.load(file1)
    data_with_site = data_with_site.fillna(0)
    mo_del.save_fingerprints_to_dataframe(data_with_site,diction_atom,128,3,True,"{}_manual_drop_duplicate_drop_atom_withtype_bond{}".format(str(128),str(3)))


    seq=sequences()
    #seq.group_seq_based_on_methylated_type()
    #seq.group_fg_based_on_methylated_type("data/input_dataframe_withoutstructure_dropatoms2048_drop_duplicate_drop_atom_withtype_bond2.csv",2048,2)
    groups=["As", "O", "S", "N", "C","Te","Se"]
    # for group in groups:
    #     parse_data.read_msa_and_encoding("{}".format(group))
    # for group in groups:
    #     build_different_input("data/group/input_dataframe_dropatoms",2048,2,seqfile="{}_onehot_encoding.csv".format(group),group=group)
    #parse_data.read_msa_and_encoding(file_name="uniprot_and_manual_align")
    # mo_del.group_by_site()
    #sepreate_input(file="data/input_data/input2048fg_dpna_bond2.csv", numbit = 2048, bond= 2)
    #create_inputdata("data/group_data",128,"6_seed")
    # data_with_site = pd.read_csv("data/seq_smiles_all_MANUAL.csv", header=0,
    #                              index_col=0)
    # site_dict = {}
    # for index in data_with_site.index:
    #     site_dict[index]=data_with_site.loc[index, "reactant_site"]
    # print(site_dict)
    # #chang the site in dictionary to manual one
    # with open("data/diction_atom_all",
    #           "wb") as dill_file:
    #     dill.dump(site_dict, dill_file)
    #mo_del.check_file_exist()
    #build_different_input(128, 2,"uniprot_and_manual_align_onehot_encoding.csv")
    # build_different_input(1024, 2)
    #parse_data.build_different_input(1024, 3)
    # print("saved input for 128bit fingerprint")
    # parse_data.build_different_input(1024,2)
    # parse_data.build_different_input(2048, 2)
    #mo_del.check_file_exist()
    #split dataset by methyl type


    # mo_del.group_by_site()
    #
    # with open('data/seq_smiles_all', 'rb') as file1:
    #     data_with_site = dill.load(file1)
    # with open('data/diction_atom_all', 'rb') as file1:
    #     diction_atom = dill.load(file1)
    #
    # # indexNames = data_with_site[data_with_site['reactant_site'] == 'NA'].index
    # # # Delete these row indexes from dataFrame
    # # data_with_site.drop(indexNames, inplace=True)
    # # #save the data after drop NA
    # # with open("data/seq_smiles_all", "wb") as dill_file:
    # #     dill.dump(data_with_site, dill_file)
    # # data_with_site.to_csv("data/diction_atom_all.csv")
    # # # print(len(data_with_site.index))
    # # # data_frame_dictionary = parse_data.group_by_domain(
    # # #     r"data\hmm_out", data_with_site)
    # # #
    # X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,128,3,
    #                             drop_atoms=True, file_name="128_drop_duplicate_withentry_drop_atom_withtype")
    # print(X)


    # X = pd.read_csv("data/input_dataframe_withoutstructure_dropatoms128_drop_duplicate_withentry_drop_atom_withtype.csv", header=0, index_col=0)
    #
    # print(X)

    'manual data'
    # with open("data/mannual_data", "rb") as dill_file:
    #     manual_data = dill.load(dill_file)
    # #save csv file to for check
    # with open("data/methyl_site_dictionary", "rb") as dill_file:
    #     methyl_site_dictionary = dill.load(dill_file)
    # #drop NA
    # indexNames = manual_data[manual_data['mainsub_mol'] == 'NA'].index
    # # Delete these row indexes from dataFrame
    # manual_data.drop(indexNames, inplace=True)
    # X = mo_del.save_fingerprints_to_dataframe(manual_data,
    #                                            methyl_site_dictionary, 2048, 3,
    #                                            drop_atoms=True,file_name="manual_2048_drop_atom_tpye")
    #

    # X = pd.read_csv("data/input_dataframe_withoutstructure_dropatomsmanual_2048_drop_atom_tpye.csv", header=0, index_col=0)
    #
    # add_dataframe = parse_data.read_msa_and_encoding("manual_align")
    #
    #
    # # add_dataframe_forPCA = add_dataframe.drop(columns=["Entry"])
    # #parse_data.create_inputdata(r"data/group_data", num_bit=128)
    #
    #
    # # # for entry in X["Entry"]:
    # # add_dataframe["Entry"]=add_dataframe["Entry"].astype(object)
    #
    # #
    # input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    # print(input_dataframe)
    # input_dataframe = input_dataframe.dropna(axis=0,how="any")
    # print(input_dataframe)
    #
    # #
    # input_dataframe.to_csv("data/input_data/input2048fg_dpna_manual.csv")
    #
    #
    



    #create_inputdata(r"data/group_data", num_bit= 128)
    #drop NA need larger memory




    # input_dataframe = pd.read_csv("data/input_data/input128fg_dpna.csv", header=0, index_col=0)
    # # with open("data/input_data/input128fg_dpna_bi_{}".format(d1), 'rb') as file1:
    # #     input_dataframe = dill.load(file1)
    # # col = [i for i in input_dataframe.columns if i not in ["Entry","label","molecular_id","methyl_type"]]
    # # input_dataframe[col] = input_dataframe[col].astype('int32')
    # # print(input_dataframe)
    # # input_dataframe = input_dataframe.reset_index()
    # # input_dataframe.drop(columns=["index"],inplace=True)
    # # print(input_dataframe)
    # # with open("data/input_data/input128fg_dpna_{}".format(d1), "wb") as dill_file:
    # #     dill.dump(input_dataframe, dill_file)
    # # input_dataframe.to_csv("data/input_data/input128fg_dpna.csv")
    # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe,i=1)
    # #protein_pca_data = copy.deepcopy(X_train).drop(columns=X_train.columns[:254])
    # #mo_del.three_D_pca(X_train,y_train,"128fg_sub_seq")

    # input_dataframe = pd.read_csv("data/input_data/input1024fg_dpna_bond2.csv", header=0, index_col=0)
    #
    #
    # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe,i=1)
    # # print(X_train)
    # # mo_del.run_PCA(X_train, y_train, "128fg_bond3")
    # # mo_del.three_D_pca(X_train, y_train, "128fg_bond3")
    # #
    # # mo_del.run_PCA(protein_pca_data, y_train, "protein_encoding_pca")
    #
    # X_train = X_train.drop(columns=["methyl_type"])
    # print(X_train.columns)
    # X_test = X_test.drop(columns=["methyl_type"])
    # y_train = y_train.drop(columns=["methyl_type"])
    # y_test = y_test.drop(columns=["methyl_type"])
    # model = mo_del.RF_model(X_train, X_test, y_train, y_test,"_input1024fg_bond2",1)
    #

if __name__ == "__main__":
    main()


