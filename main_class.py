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

    # with open('data/seq_smiles_all', 'rb') as file1:
    #     data_with_site = dill.load(file1)
    # with open('data/diction_atom_all', 'rb') as file1:
    #     diction_atom = dill.load(file1)
    #
    # indexNames = data_with_site[data_with_site['reactant_site'] == 'NA'].index
    # # Delete these row indexes from dataFrame
    # data_with_site.drop(indexNames, inplace=True)
    # print(len(data_with_site.index))
    # data_frame_dictionary = parse_data.group_by_domain(
    #     r"data\hmm_out", data_with_site)

    #X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,
    #                             drop_atoms=True, file_name="2048_drop_duplicate_withentry_drop_atom")
    #print(X)

    #X = pd.read_csv("data/input_dataframe_withoutstructure_128_drop_duplicate_withentry.csv", header=0, index_col=0)
    #X["Entry"]=X["Entry"].astype(object)

    """manual data"""
    with open("data/mannual_data", "rb") as dill_file:
        manual_data = dill.load(dill_file)
    #save csv file to for check
    with open("data/methyl_site_dictionary", "rb") as dill_file:
        methyl_site_dictionary = dill.load(dill_file)
    #drop NA
    indexNames = manual_data[manual_data['mainsub_mol'] == 'NA'].index
    # Delete these row indexes from dataFrame
    manual_data.drop(indexNames, inplace=True)
    X = mo_del.save_fingerprints_to_dataframe(manual_data,
                                               methyl_site_dictionary, 2048, 3,
                                               drop_atoms=True,file_name="manual_2048_drop_atom")

    '''
    add_dataframe = parse_data.merge_uniprot_emebeding()
    #add_dataframe = pd.read_csv("data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv", header=0, index_col=0)
    start_index=list(X.columns)[-4]
    #print(start_index)
    
    for col in add_dataframe.columns:
        if (col != "Entry") and (col != "index"):
            print(col)
            add_dataframe=add_dataframe.rename(columns={col:str(int(col)+int(start_index)+1)})
        else:
            continue
    add_dataframe.to_csv(
            "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv")
    '''
    add_dataframe = pd.read_csv(
        "data/protein_encoding/protein_seq_simple_encoding_bi_ri.csv",
        header=0, index_col=0)

    #add_dataframe = add_dataframe.rename(columns=lambda x: int(x)+int(start_index))

    # # insert_col = len(X.index)-2
    # # for entry in X["Entry"]:
    # add_dataframe["Entry"]=add_dataframe["Entry"].astype(object)
    #
    # input_dataframe_join = pd.concat([X,add_dataframe],axis=1, sort=False,ignore_index=True)
    # print(input_dataframe_join.columns.unique())
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")

    #drop NA need larger memory
    input_dataframe = input_dataframe.dropna(axis=0,how="any")
    print(input_dataframe)
    # list1=list(input_dataframe.columns)
    # list2=list(input_dataframe_join.columns)
    # for i in list1:
    #     if i not in list2:
    #         print(i)
    input_dataframe.to_csv("data/input_data/manual_2048fg_input_bi.csv")
    X_train, X_test, y_train, y_test = prepare_train_teat_data(input_dataframe)
    # model = RF_model(X_train, X_test, y_train, y_test,"_2048_drop_duplicate_double_true_label")

    # print(input_dataframe)
    #
    #
    # mo_del.run_PCA(X.iloc[:,:-1])



if __name__ == "__main__":
    main()


