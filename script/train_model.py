import dill

import parse_data
from Model_class import Model_class
import pandas as pd
import numpy as np
import copy
import time
import main_class
from datetime import date
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
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
def predict(input_dataframe):

        import joblib
        model = joblib.load('../autodata/model/rf_test_model_cvactive_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_MACCS')
        # input = copy.deepcopy(input_dataframe)
        # input =input.drop(columns=["methyl_type"])
        # input = input.drop(columns=["Entry", "molecular_id", "label"])
        input = pd.read_csv("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_MACCS_X_test.csv",header=0,index_col=0)
        input = input.drop(columns=["molecular_id","methyl_type"])
        Y = model.predict(input)
        print(Y)

        input["predict_label"] = Y
        input_frame = pd.read_csv("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_MACCS_X_test.csv",header=0,index_col=0)
        input["molecular_id"] = input_frame["molecular_id"]
        # groups=input_dataframe.groupby(["molecular_id"])
        # index = []
        # for group in groups.groups:
        #
        #         group_df=groups.get_group(group)
        #
        #         if group_df["predict_label"].sum()==1:
        #
        #                 index += list(group_df.index)
        #         else:
        #                 print(group_df)
        # input_dataframe.drop(index=index,inplace=True)
        input.to_csv("prediction_x_test.csv")
def pca_molecule(input):

        fps = input.iloc[:,:128]
        # fps["fingerprint"] = pd.DataFrame(
        #     len(fps.index) * [""])

        # for index in fps.index:
        #
        #         list1 = [str(x) for x in list(fps.loc[index])]
        #         # fps.loc[index,"fingerprint"] = "".join(list1)
        print(fps)
        # print(len(fps["fingerprint"].unique()))
        y_label = input["label"]
        fps["methyl_type"]=input["methyl_type"]
        mo_del = Model_class()
        mo_del.run_PCA(fps,y_label,"onlytmorganbit")
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

        #input_dataframe=pd.read_csv("../autodata/fingerprint/fingerprint_bit128_radius3_all_data_drop_atom_19_09.csv",header=0,index_col=0)
        #predict(input_dataframe)
        #
        # #input_dataframe=mo_del.duplicate_1_class(input_dataframe, 12)
        # #input_dataframe.drop(columns="226",inplace=True)
        # print(input_dataframe)
        # mol_id = input_dataframe["molecular_id"]
        # entry= input_dataframe["Entry"]
        # input_dataframe.drop(columns=["molecular_id","Entry"],inplace=True)
        # print(input_dataframe)
        # input_dataframe.drop_duplicates(inplace=True)
        # input_dataframe["molecular_id"]=mol_id
        # input_dataframe["Entry"] = entry
        # print(input_dataframe)
        #

        input_dataframe = pd.read_csv("../autodata/input_data/active_site/PF08241_bit_score15_coverage0.7_ACS_bit128_3_remove_redundant.csv",header=0,index_col=0)
        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
                input_dataframe)
        #Y = pd.read_csv("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_y_test.csv",header=0,index_col=0)
        #input_dataframe["label"] = Y["label"]
        X_test["label"] = y_test
        pca_molecule(X_test)

        #
        # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe)
        # #
        # #     # mo_del.three_D_pca(X_train, y_train, "{}_128_2".format(file))
        # #mo_del.run_PCA(input_dataframe, group_label, "{}".format("k_mer_sequences"))
        #
        # X_train = X_train.drop(columns=["methyl_type","molecular_id"])
        # X_test = X_test.drop(columns=["methyl_type","molecular_id"])
        # # y_train = y_train.drop(columns=["methyl_type"])
        # # y_test = y_test.drop(columns=["methyl_type"])
        # # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
        # #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
        # model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
        #                         "138fg_bond3_rf{}_{}".format(d1,"morganbit"),i=0)

if __name__ == "__main__":
    main()