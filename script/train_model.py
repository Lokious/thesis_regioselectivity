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
from molecular_class import Molecule
from rdkit import Chem
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
                       "{}".format("128bitfp"))
def predict(name="N_AA_properties_encoding_MACCSkey"):

        import joblib
        print(name)
        model = joblib.load('../autodata/model/rf_test_model_cv{}'.format(name))
        # input = copy.deepcopy(input_dataframe)
        # input =input.drop(columns=["methyl_type"])
        # input = input.drop(columns=["Entry", "molecular_id", "label"])
        input = pd.read_csv("../autodata/model/{}_X_test.csv".format(name),header=0,index_col=0)
        input = input.drop(columns=["molecular_id","methyl_type","atom_index"])
        Y = model.predict(input)
        #print(Y)
        y_test= pd.read_csv("../autodata/model/{}_y_test.csv".format(name),header=0,index_col=0)
        input["predict_label"] = Y
        input["true_label"]=y_test["label"]
        input_frame = pd.read_csv("../autodata/model/{}_X_test.csv".format(name),header=0,index_col=0)
        input["molecular_id"] = input_frame["molecular_id"]
        input["atom_index"] = input_frame["atom_index"]
        input["methyl_type"] = input_frame["methyl_type"]
        #print(input["methyl_type"].unique())
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
        input.to_csv("{}prediction_x_test.csv".format(name))
def pca_molecule(input):

        fps = input.iloc[:,:4096]
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
        mo_del.run_PCA(fps,y_label,"2048molecules")
def rename_column():
        input_dataframe = pd.read_csv(
                "../autodata/fingerprint/fingerprint_bit128_radius3_all_data.csv",
                header=0, index_col=0)
        # other_info = pd.DataFrame(input_dataframe,columns=['Entry','atom_index','label','methyl_type','molecular_id'])

        print(input_dataframe)
        molecule = Molecule()

        map_dict = {}
        list_columns = []
        for key in range(128):
                # print(key)
                map_dict[str(key)] = "molecule_" + str(key)
                list_columns.append("molecule_" + str(key))

        for i in list(range(128, 256)):
                map_dict[str(i)] = "atom_" + str(i)
                list_columns.append("atom_" + str(i))
        print(map_dict)
        fingerprint_df = input_dataframe.rename(columns=map_dict)
        fingerprint_df.dropna(inplace=True)
        print(fingerprint_df)
        fingerprint_df.to_csv(
                "../autodata/fingerprint/fingerprint_bit128_radius3_all_data.csv")
        columns = list_columns
        input_dataframe = fingerprint_df.drop_duplicates(subset=columns)
        input_dataframe = input_dataframe.dropna()
        print(input_dataframe)
        input_dataframe.to_csv(
                "../autodata/input_data/fingerprint_bit128_radius3_all_data_morganfingerprint.csv")
        #rename for maccs key
        input_dataframe = pd.read_csv(
                "../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data_2_11.csv",
                header=0, index_col=0)
        # other_info = pd.DataFrame(input_dataframe,columns=['Entry','atom_index','label','methyl_type','molecular_id'])
        input_dataframe = input_dataframe.drop(columns=["0", "167"])
        print(input_dataframe)
        molecule = Molecule()
        maccs_dict = molecule.MACCSkey_dictionary
        map_dict = {}
        for key in maccs_dict.keys():
                # print(key)
                map_dict[str(key)] = maccs_dict[key][0]

        for i in list(range(168, 334)):
                map_dict[str(i)] = "atom" + maccs_dict[(i - 167)][0]
        print(map_dict)
        fingerprint_df = input_dataframe.rename(columns=map_dict)
        print(input_dataframe)
        fingerprint_df.dropna(inplace=True)
        print(fingerprint_df)
        fingerprint_df.to_csv(
                "../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data.csv")

        input_dataframe = pd.read_csv(
                "../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data_2_11.csv",
                header=0, index_col=0)
        # other_info = pd.DataFrame(input_dataframe,columns=['Entry','atom_index','label','methyl_type','molecular_id'])

        print(input_dataframe)
        columns = [str(x) for x in list(range(1, 334))]
        input_dataframe = input_dataframe.drop_duplicates(subset=columns)
        input_dataframe = input_dataframe.dropna()
        print(input_dataframe)
        input_dataframe = input_dataframe.drop(columns=["0", "167"])
        input_dataframe = input_dataframe.rename(columns=map_dict)
        input_dataframe.dropna(inplace=True)
        print(input_dataframe)
        input_dataframe.to_csv(
                "../autodata/input_data/MACCS_fingerprint_bit167_radius3_all_data_2_11_change_name.csv")
        # input_dataframe=pd.read_csv("../autodata/input_data/fingerprint_bit128_radius3_all_data_morganfingerprint.csv")
        print(input_dataframe)

def use_manual_data_for_test(data="../autodata/mannual_data.csv"):
        dataset = pd.read_csv(data,index_col=0,header=0)
        dataset['substrate_mol']= pd.DataFrame(len(dataset.index) * [0]).astype('object')
        dataset=dataset.drop(columns=["PDB"])
        dataset.dropna(inplace=True)
        for i in dataset.index:
                sub_smile = dataset.loc[i,"main_sub"]
                mol = Chem.MolFromSmiles(sub_smile)
                dataset.loc[i,'substrate_mol'] = mol
        mo_del = Model_class()
        fg=mo_del.MACCSkey_fingerprint_df_preapare(substrate_mol_df=dataset,radius = 3)
        print(fg)
        import joblib
        name="MACCS_fingerprint_bit167_radius3_all_data_no_same_sub_1"
        model = joblib.load('../autodata/model/rf_test_model_cv{}'.format(name))
        # input = copy.deepcopy(input_dataframe)
        # input =input.drop(columns=["methyl_type"])
        # input = input.drop(columns=["Entry", "molecular_id", "label"])
        input =copy.deepcopy(fg).drop(columns=[0,167])
        print(input["label"].sum())
        input = input.drop(columns=["molecular_id","atom_index","label"])
        Y = model.predict(input)
        print(Y)
        fg['predict']= pd.DataFrame(len(input.index) * [0]).astype('string')
        fg['predict']=Y
        count =0
        for index in fg.index:
                if fg.loc[index,'predict']== fg.loc[index,"label"]:
                        count +=1
        print(count)
        print((count/len(input.index))) #95.37% accuracy probabily all the substrate has already present in trainning data

def main():

        mo_del = Model_class()
        today = date.today()
        # dd/mm/YY
        d1 = today.strftime("%d_%m_%Y")
        #
        # input_dataframe = pd.read_csv(
        #         "../autodata/fingerprint/fingerprint_bit2048_radius3_all_data.csv",
        #         header=0, index_col=0)
        # # other_info = pd.DataFrame(input_dataframe,columns=['Entry','atom_index','label','methyl_type','molecular_id'])
        #
        # print(input_dataframe)
        # molecule = Molecule()
        #
        # map_dict = {}
        # list_columns = []
        # for key in range(2048):
        #         # print(key)
        #         map_dict[str(key)] = "molecule_" + str(key)
        #         list_columns.append("molecule_" + str(key))
        #
        # for i in list(range(2048,4096)):
        #         map_dict[str(i)] = "atom_" + str(i)
        #         list_columns.append("atom_" + str(i))
        # print(map_dict)
        # fingerprint_df = input_dataframe.rename(columns=map_dict)
        # fingerprint_df.dropna(inplace=True)
        # print(fingerprint_df)
        # fingerprint_df.to_csv(
        #         "../autodata/fingerprint/fingerprint_bit2048_radius3_all_data.csv")
        # df_1024 = pd.read_csv(
        #         "../autodata/fingerprint/fingerprint_bit2048_radius3_all_data.csv",
        #         header=0, index_col=0)
        # df_1024.dropna(inplace=True)
        # pca_molecule(input=df_1024)
        #use_manual_data_for_test(data="../autodata/mannual_data.csv")
        # for coverage in [50, 60, 70, 80, 90]:
        #         for bitscore in [11, 15]:
        #                 predict(name=
        #                 "active_site_128fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_{}_{}_MACCS_no_same_sub".format(
        #                 bitscore, coverage))
        #predict(name="6_seed_onehot_encoding_MACCSkey")
        # input = pd.read_csv("../autodata/input_data/active_site/C_seed_onehot_encoding_sepreate_MACCSkey.csv",header=0,index_col=0)
        # pca_molecule(input)
        # input = pd.read_csv("../autodata/fingerprint/fingerprint_bit128_radius3_all_data.csv",header=0,index_col=0)
        # pca_molecule(input)
        # predict(name="N_AA_properties_encoding_MACCSkey_no_same_sub")
        # parse_data.molecular_accuracy("N_AA_properties_encoding_MACCSkey_no_same_subprediction_x_test.csv")
        # predict(name="O_AA_properties_encoding_MACCSkey_no_same_sub")
        # parse_data.molecular_accuracy(
        #         "O_AA_properties_encoding_MACCSkey_no_same_subprediction_x_test.csv")
        # predict(name="N_AA_properties_encoding_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "N_AA_properties_encoding_MACCSkeyprediction_x_test.csv")
        # predict(name="O_AA_properties_encoding_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "O_AA_properties_encoding_MACCSkeyprediction_x_test.csv")
        # predict(name="onlyfingerprint_bit128_radius3_all_data_morganfingerprint_no_same_sub")
        # parse_data.molecular_accuracy(
        #         "onlyfingerprint_bit128_radius3_all_data_morganfingerprint_no_same_subprediction_x_test.csv")
        # predict(name="onlyfingerprint_bit128_radius3_all_data_morganfingerprint")
        # parse_data.molecular_accuracy(
        #         "onlyfingerprint_bit128_radius3_all_data_morganfingerprintprediction_x_test.csv")
        # predict(name="MACCS_fingerprint_bit167_radius3_all_data_2_11_change_name")
        # parse_data.molecular_accuracy(
        #         "MACCS_fingerprint_bit167_radius3_all_data_2_11_change_nameprediction_x_test.csv")

        # #rename_column()
        #
        # # N methyltransferase
        # X = pd.read_csv(
        #  "../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data.csv",
        #  header=0, index_col=0)
        # add_dataframe = pd.read_csv(
        #  "../autodata/protein_encoding/k_mer/C_k_mer_encoding_without_align_26_08.csv",
        #  header=0, index_col=0)
        # #add_dataframe["Entry"] = add_dataframe.index
        # #add_dataframe.reset_index(drop=True, inplace=True)
        # print(add_dataframe)
        # # drop all zero columns
        # add_dataframe = add_dataframe.drop(columns=["group","ID"])
        # add_dataframe = (add_dataframe.loc[:, add_dataframe.sum() != 0.0])
        # add_dataframe.drop_duplicates(subset="Entry", inplace=True)
        # input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
        # print(input_dataframe)
        # input_dataframe = input_dataframe.dropna(axis=0, how="any")
        # # input_dataframe["methyl_type"] = input_dataframe["group"]
        #
        # columns=input_dataframe.select_dtypes(include=['float','int64']).columns
        # input_dataframe[columns]=input_dataframe[columns].astype('int32')
        # print(input_dataframe)
        # input_dataframe.to_csv(
        #  "../autodata/input_data/16C_k_mer_encoding_without_align_MACCSkey.csv")
        input_dataframe = pd.read_csv("../autodata/input_data/active_site/PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS.csv",header=0,index_col=0)
        print(input_dataframe)
        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
                input_dataframe)
        X_train = X_train.drop(
                columns=["methyl_type", "molecular_id", "atom_index"])
        # save x test for further analysis result
        y_test.to_csv(
                "../autodata/model/PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS_y_test.csv")
        X_test.to_csv(
                "../autodata/model/PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS_X_test.csv")
        X_test = X_test.drop(
                columns=["methyl_type", "molecular_id", "atom_index"])
        model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                                 "PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS", i=0)
        parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
                              input_data="../autodata/input_data/active_site/PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS.csv")
        predict(name="PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCS")
        parse_data.molecular_accuracy(
                "PF08241_bit_score11_coverage0.6_ACS_bit167_3_remove_redundant_MACCSprediction_x_test.csv")
        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        # #                       input_data="../autodata/input_data/16C_k_mer_encoding_without_align_MACCSkey.csv")
        # predict(name="16C_k_mer_encoding_without_align_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "16C_k_mer_encoding_without_align_MACCSkeyprediction_x_test.csv")
        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="../autodata/input_data/N_k_mer_encoding_without_align_MACCSkey.csv")
        # predict(name="N_k_mer_encoding_without_align_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "N_k_mer_encoding_without_align_MACCSkeyprediction_x_test.csv")
        # parse_data.read_msa_and_encoding("6_seed")

        #predict("")
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

        #input_dataframe=pd.read_csv("../autodata/input_data/MACCS_fingerprint_bit167_radius3_all_data_31_10.csv",header=0,index_col=0)
        # input_dataframe=pd.read_csv(r"E:\Download\regioselectivity_prediction\autodata\input_data\active_site\O_AA_properties_encoding_MACCSkey.csv",header=0,index_col=0)
        # label=pd.read_csv(r"E:\Download\regioselectivity_prediction\autodata\input_data\active_site\O_AA_properties_encoding_MACCSkey.csv",header=0,index_col=0)
        # input_dataframe["label"]=label["label"]
        # print(input_dataframe)
        # pca_molecule(input_dataframe)


        # predict("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_MACCS")
        # predict("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_50_MACCS")
        # predict("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_11_70_MACCS")
        # predict("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_11_50_MACCS")


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
        '''
        input_dataframe = pd.read_csv("../autodata/input_data/active_site/PF08241_bit_score15_coverage0.7_ACS_bit128_3_remove_redundant.csv",header=0,index_col=0)
        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
                input_dataframe)
        #Y = pd.read_csv("active_site_167fg_bi_type_bond3_rf_PF08241_ACS_remove_redundant_15_70_y_test.csv",header=0,index_col=0)
        #input_dataframe["label"] = Y["label"]
        X_test["label"] = y_test
        pca_molecule(X_test)
        '''
        #######train on only morgan fp#########
        # input_dataframe=pd.read_csv( "../autodata/input_data/fingerprint_bit128_radius3_all_data_morganfingerprint.csv",header=0,index_col=0)
        # print(input_dataframe)
        # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        #         input_dataframe,group_column="main_sub")
        # y_test.to_csv(
        #         "../autodata/model/onlyfingerprint_bit128_radius3_all_data_morganfingerprint_no_same_sub_y_test.csv")
        # X_test.to_csv(
        #         "../autodata/model/onlyfingerprint_bit128_radius3_all_data_morganfingerprint_no_same_sub_X_test.csv")
        # X_train = X_train.drop(
        #         columns=["methyl_type", "molecular_id", "atom_index"])
        # # save x test for further analysis result
        #
        # X_test = X_test.drop(
        #         columns=["methyl_type", "molecular_id", "atom_index"])
        # model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
        #                          "onlyfingerprint_bit128_radius3_all_data_morganfingerprint_no_same_sub", i=0)

        # input_dataframe=pd.read_csv( "../autodata/input_data/MACCS_fingerprint_bit167_radius3_all_data_2_11.csv",header=0,index_col=0)
        # print(input_dataframe)
        # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        #         input_dataframe,group_column="main_sub",i=1)
        # print(X_test)
        #
        # y_test.to_csv(
        #         "../autodata/model/MACCS_fingerprint_bit167_radius3_all_data_no_same_sub_1_y_test.csv")
        # X_test.to_csv(
        #         "../autodata/model/MACCS_fingerprint_bit167_radius3_all_data_no_same_sub_1_X_test.csv")
        # X_train = X_train.drop(
        #         columns=["methyl_type", "molecular_id", "atom_index"])
        # # save x test for further analysis result
        #
        # X_test = X_test.drop(
        #         columns=["methyl_type", "molecular_id", "atom_index"])
        # model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
        #                          "MACCS_fingerprint_bit167_radius3_all_data_no_same_sub_1", i=0)
        # #
        # predict(name="O_seed_onehot_encoding_sepreate_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "O_seed_onehot_encoding_sepreate_MACCSkeyprediction_x_test.csv")
        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="../autodata/input_data/all_k_mer_encoding_MACCSkey.csv")
        # predict(name="O_seed_onehot_encoding_sepreate_MACCSkey_no_same_sub")
        # parse_data.molecular_accuracy(
        #         "O_seed_onehot_encoding_sepreate_MACCSkey_no_same_subprediction_x_test.csv")
        parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
                              input_data="../autodata/input_data/active_site/PF08241_bit_score11_coverage0.6_ACS_bit128_3_remove_redundant_MACCS.csv")
        predict(name="PF08241_bit_score11_coverage0.6_ACS_bit128_3_remove_redundant_MACCS")
        parse_data.molecular_accuracy(
                "PF08241_bit_score11_coverage0.6_ACS_bit128_3_remove_redundant_MACCSprediction_x_test.csv")
        # predict(name="S_AA_properties_encoding_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "S_AA_properties_encoding_MACCSkeyprediction_x_test.csv")
        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="../autodata/input_data/O_k_mer_encoding_without_align_MACCSkey.csv")

        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="../autodata/input_data/C_k_mer_encoding_without_align_MACCSkey.csv")
        # predict(name="S_k_mer_encoding_without_align_MACCSkey")
        # parse_data.molecular_accuracy(
        #         "S_k_mer_encoding_without_align_MACCSkeyprediction_x_test.csv")
        # parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="../autodata/input_data/S_k_mer_encoding_without_align_MACCSkey.csv")
        # predict(name="O_AA_encoding_without_align_MACCSkey_fromstructure")
        # parse_data.molecular_accuracy(
        #         "O_AA_encoding_without_align_MACCSkey_fromstructureprediction_x_test.csv")
        # train_sub =parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="traindataO_AA_MACCS.csv")
        # test_sub =parse_data.molecule_number_count(substrate_df="../autodata/seq_smiles_all.csv",
        #                       input_data="testdataO_AA_MACCS.csv")
        # subs = set(list(test_sub) +list(train_sub))
        # print(subs)
        # print(len(subs))
        # test = pd.read_csv("testdataN_AA.csv",header=0,index_col=0)
        # print(test)
        # test = test.dropna()
        # print(test)
        # X_test = (copy.deepcopy(test)).drop(
        #         columns=["label", "invariants_kmer", "Entry",
        #                  "invariants_radius",
        #                   "structure"])
        # X_test.to_csv("../autodata/model/{}_X_test.csv".format("166fg_rf11_6_N_AA_with_structure"))
        # y_test = test["label"]
        # y_test.to_csv("../autodata/model/{}_y_test.csv".format(
        #         "166fg_rf11_6_N_AA_with_structure"))
        #
        # predict(name="166fg_rf11_6_N_AA_with_structure")
        # parse_data.molecular_accuracy(
        #         "166fg_rf11_6_N_AA_with_structureprediction_x_test.csv")

        #
        #
        # predict(name="6_seed_onehot_encoding_MACCSkey_no_same_sub")
        # parse_data.molecular_accuracy(
        #         "6_seed_onehot_encoding_MACCSkey_no_same_subprediction_x_test.csv")
        # parse_data.molecule_number_count(
        #         substrate_df="../autodata/seq_smiles_all.csv",
        #         input_data="../autodata/input_data/active_site/6_seed_onehot_encoding_MACCSkey.csv")
        # X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(input_dataframe)
        # #
        # #     # mo_del.three_D_pca(X_train, y_train, "{}_128_2".format(file))
        # #mo_del.run_PCA(input_dataframe, group_label, "{}".format("k_mer_sequences"))
        # X_test.to_csv("X_test_for_only_maccskeys.csv")
        # y_test.to_csv("y_test_for_only_maccskeys.csv")
        # X_train = X_train.drop(columns=["methyl_type","molecular_id","atom_index"])
        # X_test = X_test.drop(columns=["methyl_type","molecular_id","atom_index"])
        # print(X_test)
        # # y_train = y_train.drop(columns=["methyl_type"])
        # # y_test = y_test.drop(columns=["methyl_type"])
        # # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
        # #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
        # model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
        #                         "166fg_bond3_rf{}_{}".format(d1,"MACCS_keys"),i=0)

if __name__ == "__main__":
    main()