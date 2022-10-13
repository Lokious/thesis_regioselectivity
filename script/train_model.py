import dill

import parse_data
from Model_class import Model_class
import pandas as pd
import numpy as np
import copy
import time
import main_class
from datetime import date
def run_pca_for_kmer_sequences():
        mo_del = Model_class()
        sequence_data = pd.read_csv(
                "../autodata/protein_encoding/k_mer/all_k_mer_encoding_sepreate_without_align.csv",
                header=0, index_col=0)
        sequence_data['methyl_type'] = sequence_data["group"]
        group_label = sequence_data['methyl_type']
        sequence_data.drop("Entry", axis=1, inplace=True)
        sequence_data.drop("group", axis=1, inplace=True)
        mo_del.run_PCA(sequence_data, group_label,
                       "{}".format("k_mer_sequences"))

def main():

        mo_del = Model_class()
        today = date.today()
        # dd/mm/YY
        d1 = today.strftime("%d_%m_%Y")

        #parse_data.read_msa_and_encoding("6_seed")


        #sequence_data = pd.read_csv("../autodata/input_data/input128fg_dpna_bond3_O_seed_onehot_encoding.csv.csv",header=0,index_col=0)
        # sequence_data['methyl_type']=sequence_data["group"]
        # # print(sequence_data.iloc[:,256:])
        # # sequence_data=sequence_data.iloc[:,256:]
        # group_label=sequence_data['methyl_type']
        # #sequence_data.drop("methyl_type",axis=1,inplace=True)
        # sequence_data.drop("Entry",axis=1,inplace=True)
        # sequence_data.drop("group",axis=1,inplace=True)
        #sequence_data.drop("molecular_id",axis=1,inplace=True)
        #print(sequence_data)

        #mo_del.hierarchical_clustering(sequence_data,group_label)
        # filename_list=[ "S","O","N", "C"]
        # for file in filename_list:
        #     input_dataframe = pd.read_csv("../autodata/input_data/bit_info/input128fg_dpna_bond3_{}_seed_onehot_encoding.csv.csv".format(file), header=0, index_col=0)
        #     #input_dataframe.dropna(inplace=True)

        input_dataframe=pd.read_csv("../autodata/input_data/active_site/PF08241PF01795_bit_score15_coverage0.8_ACS_bit128_3_remove_redundant.csv",header=0,index_col=0)
        #input_dataframe=mo_del.duplicate_1_class(input_dataframe, 12)
        #input_dataframe.drop(columns="226",inplace=True)
        print(input_dataframe)
        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe)
        #
        #     # mo_del.three_D_pca(X_train, y_train, "{}_128_2".format(file))
        #mo_del.run_PCA(sequence_data, group_label, "{}".format("k_mer_sequences"))

        X_train = X_train.drop(columns=["methyl_type"])
        X_test = X_test.drop(columns=["methyl_type"])
        y_train = y_train.drop(columns=["methyl_type"])
        y_test = y_test.drop(columns=["methyl_type"])
        # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
        #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
        model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                                "active_site_128fg_bi_type_bond3_rf{}_{}".format(d1,"PF08241PF01795_bit_score15_coverage80"),i=0)

if __name__ == "__main__":
    main()