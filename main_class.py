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
def main():

    mo_del = Model_class()
    #mo_del.check_file_exist()
    #split dataset by methyl type
    #mo_del.group_by_site()
    #
    with open('data/seq_smiles_all', 'rb') as file1:
        data_with_site = dill.load(file1)
    with open('data/diction_atom_all', 'rb') as file1:
        diction_atom = dill.load(file1)

    indexNames = data_with_site[data_with_site['reactant_site'] == 'NA'].index
    # Delete these row indexes from dataFrame
    data_with_site.drop(indexNames, inplace=True)
    #save the data after drop NA
    with open("data/seq_smiles_all", "wb") as dill_file:
        dill.dump(data_with_site, dill_file)
    data_with_site.to_csv("data/diction_atom_all.csv")
    # # print(len(data_with_site.index))
    # # data_frame_dictionary = parse_data.group_by_domain(
    # #     r"data\hmm_out", data_with_site)
    # #
    X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,512,3,
                                drop_atoms=True, file_name="512_drop_duplicate_withentry_drop_atom_withtype")
    print(X)

    X = pd.read_csv("data/input_dataframe_withoutstructure_512_drop_duplicate_withentry_drop_atom.csv", header=0, index_col=0)


    # """manual data"""
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
    #
    #X = pd.read_csv("data/input_dataframe_withoutstructure_dropatomsmanual_2048_drop_atom_tpye.csv", header=0, index_col=0)
    #
    # add_dataframe = parse_data.merge_uniprot_emebeding()
    # add_dataframe = pd.read_csv("data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv", header=0, index_col=0)
    # print(add_dataframe)
    # start_index=list(X.columns)[-5]
    # #print(start_index)

    # for col in add_dataframe.columns:
    #     if (col != "Entry") and (col != "index"):
    #         print(col)
    #         add_dataframe=add_dataframe.rename(columns={col:str(int(col)+int(start_index)+1)})
    #     else:
    #         continue
    # add_dataframe.to_csv(
    #         "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv")

    # add_dataframe = pd.read_csv(
    #     "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv",
    #     header=0, index_col=0)
    # add_dataframe_forPCA = add_dataframe.drop(columns=["Entry"])



    # # for entry in X["Entry"]:
    # add_dataframe["Entry"]=add_dataframe["Entry"].astype(object)
    #
    # input_dataframe_join = pd.concat([X,add_dataframe],axis=1, sort=False,ignore_index=True)
    # print(input_dataframe_join.columns.unique())
    #
    # input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    #
    # #drop NA need larger memory
    # input_dataframe = input_dataframe.dropna(axis=0,how="any")
    # print(input_dataframe)
    # #
    # input_dataframe.to_csv("data/input_data/input128fg_dpna_bi.csv")
    #parse_data.create_inputdata(r"data/group_data", num_bit= 2048)








    '''
    input_dataframe = pd.read_csv("data/input_data/manual_2048fg_input_bi_type.csv", header=0, index_col=0)

    # col = [i for i in input_dataframe.columns if i not in ["Entry","label","molecular_id","methyl_type"]]
    #
    # input_dataframe[col] = input_dataframe[col].astype('int32')
    # input_dataframe.to_csv("data/input_data/input128fg_dpna_bi.csv")
    # print(input_dataframe)
    #input_dataframe = input_dataframe.reset_index()

    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe)

    mo_del.run_PCA(X_train)
    X_train = X_train.drop(columns=["methyl_type"])
    X_test = X_test.drop(columns=["methyl_type"])
    y_train = y_train.drop(columns=["methyl_type"])
    y_test = y_test.drop(columns=["methyl_type"])
    model = mo_del.RF_model(X_train, X_test, y_train, y_test,"_manual_2048fg_input_bi_type")
    '''
    # print(input_dataframe)
    #
    #




if __name__ == "__main__":
    main()


