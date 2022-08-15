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
            input_data = pd.read_csv(
                "data/input_data/input{0}fg_dpna_{1}.csv".format(
                    str(num_bit),
                    name1), header=0, index_col=0)
            print(input_data)
        except:
            print("data/input_data/input{0}fg_dpna_{1}.csv is missing".format(
                    str(num_bit),
                    name1))
            try:
                print("loading data:{}-------".format(name1))
                X = pd.read_csv(
                    "data/input_dataframe_withoutstructure_dropatoms{}_{}_withentry_drop_atom_type.csv".format(
                        str(num_bit),
                        name1), header=0, index_col=0)
                print(X)
            except:
                print("data/input_dataframe_withoutstructure_dropatoms{}_{}_withentry_drop_atom_type.csv is missing".format(
                        str(num_bit),
                        name1))
                print("creat fingerprint {}-------".format(name1))
                with open('{}/seq_smiles_all[{}'.format(directory, name1),
                          'rb') as file1:
                    data_with_site = dill.load(file1)
                with open('{}/diction_atom_all[{}'.format(directory, name1),
                          'rb') as file2:
                    diction_atom = dill.load(file2)
                print(diction_atom)
                X = mo_del.save_fingerprints_to_dataframe(data_with_site,
                                                          diction_atom,
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
                                       num_bit=num_bit, add_dataframe="data/protein_encoding/6_seed_onehot_encoding_rn_{}fg.csv".format(str(num_bit)))
            input_data = pd.read_csv("data/input_data/input{0}fg_dpna_{1}.csv".format(str(num_bit),
                                                            name1), header=0, index_col=0)
        print(input_data)
        if len(input_data.index)>0:
            run_model_for_group_data(
                input_data, filename=name1, num_bit=num_bit)
        else:
            print("Empty dataframe: {}".format(name1))
def merge_sequence_and_fingerprint(x:pd.DataFrame,filename:str="",num_bit:int=2048, add_dataframe=""):

    if add_dataframe == "":
        add_dataframe = pd.read_csv(
            "data/protein_encoding/6_seed_onehot_encoding.csv",
            header=0, index_col=0)
    else:
        try:
            add_dataframe = pd.read_csv("{}".format(add_dataframe), header=0, index_col=0)
        except:
            add_dataframe = pd.read_csv(
                "data/protein_encoding/6_seed_onehot_encoding.csv",
                header=0, index_col=0)

    print(add_dataframe)
    start_index = (num_bit*2-1)
    #print(start_index)
    if list(add_dataframe.columns)[0] != str(int(start_index)+1):
        map_dictionary = {}
        for col in add_dataframe.columns:
            if (col != "Entry") and (col != "index"):

                map_dictionary[col] = str(int(col) + int(start_index))

            else:
                continue
        add_dataframe = add_dataframe.rename(columns=map_dictionary)
    add_dataframe.to_csv(
            "data/protein_encoding/6_seed_onehot_encoding_rn_{}fg.csv".format(str(num_bit)))
    add_dataframe = pd.read_csv(
        "data/protein_encoding/6_seed_onehot_encoding_rn_{}fg.csv".format(str(num_bit)),
        header=0, index_col=0)
    input_dataframe = x.merge(add_dataframe, on="Entry", how="left")

    #drop NA need larger memory
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    col = [i for i in input_dataframe.columns if
           i not in ["Entry", "label", "molecular_id", "methyl_type"]]
    input_dataframe[col] = input_dataframe[col].astype('int32')
    input_dataframe = (input_dataframe.reset_index()).drop(columns=["index"])
    print(input_dataframe)
    #
    input_dataframe.to_csv("data/input_data/input{0}fg_dpna_{1}.csv".format(str(num_bit),filename))
    with open("data/input_data/input{0}fg_dpna_{1}".format(str(num_bit),filename), "wb") as dill_file:
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

def build_different_input(num_bit:int=0,radius:int=0):
    """

    :param num_bit:
    :param radius:
    :return:
    """
    mo_del = Model_class()

    try:
        X = pd.read_csv(
            "data/input_dataframe_withoutstructure_dropatoms{}_drop_duplicate_drop_atom_withtype_bond{}.csv".format(
                str(num_bit), str(radius)), header=0, index_col=0)
    except:
        data_with_site = pd.read_csv("data/seq_smiles_all_MANUAL.csv", header=0, index_col=0)
        with open('data/diction_atom_all', 'rb') as file1:
            diction_atom = dill.load(file1)
        X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,
                                                  num_bit, radius,
                                                  drop_atoms=True,
                                                  file_name="{}_drop_duplicate_drop_atom_withtype_bond{}".format(str(num_bit),str(radius)))
    print(X)
    add_dataframe = pd.read_csv("data/protein_encoding/6_seed_onehot_encoding.csv", header=0, index_col=0)
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
                "data/protein_encoding/6_seed_onehot_encoding_rn_{}fg.csv".format(str(num_bit)))
    print("merging fingerprint and sequences encoding---------")
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    print(input_dataframe)
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    print(input_dataframe)
    print("saving input data.......")
    #
    input_dataframe.to_csv("data/input_data/input{}fg_dpna_bond{}.csv".format(str(num_bit),str(radius)))


def main():
    today = date.today()
    # dd/mm/YY
    d1 = today.strftime("%d_%m_%Y")
    mo_del = Model_class()
    # mo_del.group_by_site()
    #create_inputdata("data/group_data",512)
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
    #build_different_input(512, 2)
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
    



    create_inputdata(r"data/group_data", num_bit= 128)
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

    input_dataframe = pd.read_csv("data/input_dataframe_withoutstructure_dropatoms128_drop_duplicate_drop_atom_withtype_bond3.csv", header=0, index_col=0)
    #
    #
    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe,i=1)
    # print(X_train)
    # mo_del.run_PCA(X_train, y_train, "128fg_bond3")
    # mo_del.three_D_pca(X_train, y_train, "128fg_bond3")

    #mo_del.run_PCA(protein_pca_data, y_train, "protein_encoding_pca")
    #
    # X_train = X_train.drop(columns=["methyl_type"])
    # print(X_train.columns)
    # X_test = X_test.drop(columns=["methyl_type"])
    # y_train = y_train.drop(columns=["methyl_type"])
    # y_test = y_test.drop(columns=["methyl_type"])
    model = mo_del.RF_model(X_train, X_test, y_train, y_test,"_input128fg_type_withoutstructure{}".format(1),1)


if __name__ == "__main__":
    main()


