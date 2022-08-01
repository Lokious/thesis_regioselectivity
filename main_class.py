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

def main():

    mo_del = Model_class()
    '''
    with open('data/seq_smiles_all', 'rb') as file1:
        data_with_site = dill.load(file1)
    with open('data/diction_atom_all', 'rb') as file1:
        diction_atom = dill.load(file1)

    indexNames = data_with_site[data_with_site['reactant_site'] == 'NA'].index
    # Delete these row indexes from dataFrame
    data_with_site.drop(indexNames, inplace=True)
    # print(len(data_with_site.index))
    # data_frame_dictionary = parse_data.group_by_domain(
    #     r"data\hmm_out", data_with_site)
    
    X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom, 128,
                                   3, file_name="128_drop_duplicate_withentry")
    '''
    X = pd.read_csv("data/input_dataframe_withoutstructure_128_drop_duplicate_withentry.csv", header=0, index_col=0)
    print(X)
    add_dataframe = parse_data.read_msa_and_encoding()
    input_dataframe = pd.concat([X, add_dataframe], axis=
    0)
    print(input_dataframe)


    mo_del.run_PCA(X.iloc[:,:-1])



if __name__ == "__main__":
    main()


