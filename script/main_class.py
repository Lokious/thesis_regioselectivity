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

import parse_data
from Model_class import Model_class
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
from sequence import Sequences


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

def build_different_input(auto="",x="",num_bit:int=0,radius:int=0,seqfile:str="6_seed_onehot_encoding.csv",group=""):
    """

    :param num_bit:
    :param radius:
    :return:
    """
    mo_del = Model_class()
    try:
        input_dataframe=pd.read_csv(
            "../{}data/input_data/bit_info/input{}fg_dpna_bond{}_{}.csv".format(auto,
                                                                       str(num_bit),
                                                                       str(radius),
                                                                       seqfile),header=0,index_col=0)
        print(input_dataframe)
    except:
        print("creating inputframe for :{}".format(seqfile))
        try:
            if x=="":
                X = pd.read_csv(
                    "../{}data/input_dataframe_withoutstructure_dropatoms{}_drop_duplicate_drop_atom_withtype_bond{}.csv".format(auto,
                        str(num_bit), str(radius)), header=0, index_col=0)
            else:
                X = pd.read_csv("{}".format(x),header=0,index_col=0)
        except:
            if x=="":
                print("meet problem while trying to reading X")
                data_with_site = pd.read_csv("../{}data/seq_smiles_all.csv".format(auto), header=0, index_col=0)
                with open('../{}data/diction_atom_all'.format(auto), 'rb') as file1:
                    diction_atom = dill.load(file1)
                create_fingerprint=input("please input Y to continue:")
                if create_fingerprint=="Y" or "y":
                    X = mo_del.save_fingerprints_to_dataframe(data_with_site, diction_atom,
                                                          num_bit, radius,
                                                          drop_atoms=True,
                                                          file_name="drop_duplicate_drop_atom_withtype".format(str(num_bit),str(radius)))
                else:
                    print("existing-----")
                    exit()
            else:
                raise IOError("The input x is not exit")

        print(X)
        X.dropna(inplace=True)
        print("x after drop na:")
        print(X)
        try:
            add_dataframe = pd.read_csv("../{}data/protein_encoding/{}_{}fg_rm.csv".format(auto,seqfile,str(num_bit)), header=0, index_col=0)
            print(add_dataframe)
        except:
            print("../{}data/protein_encoding/{}_{}fg.csv missing, build protein_encoding data------".format(auto,seqfile,str(num_bit)))
            # create_add_dataframe=input("please input Y to continue:")
            # if create_add_dataframe == "y" or "Y":
            add_dataframe=parse_data.read_msa_and_encoding(group)
            # else:
            #     print("existing-----")
            #     exit()
        start_index = num_bit*2
        #print(start_index)
        if list(add_dataframe.columns)[0] != str(start_index):
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
                    "../{}data/protein_encoding/{}_{}fg_rm.csv".format(auto,seqfile,str(num_bit)))
        print("merging fingerprint and sequences encoding---------")
        input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
        print(input_dataframe)
        input_dataframe = input_dataframe.dropna(axis=0, how="any")
        print(input_dataframe)
        print("saving input data.......")
        #
        input_dataframe.to_csv("../{}data/input_data/bit_info/input{}fg_dpna_bond{}_{}.csv".format(auto,str(num_bit),str(radius),seqfile))
def use_k_merencoding_for_create_input(x="",num_bit:int=0,radius:int=0,seqfile:str="all_k_mer_encoding_sepreate_without_align.csv",group=""):
    """
    This function is to combine k_mer protein encoding and substrate fingerprint

    :param x:
    :param num_bit:
    :param radius:
    :param seqfile:
    :param group:
    :return:
    """
    print(x)
    #covert dtype to sting, avoid incorrect NA while merging
    X = pd.read_csv("{}".format(x), header=0, index_col=0)
    X["Entry"]=(X["Entry"].astype("string"))
    print(X.dtypes)
    add_dataframe=pd.read_csv("../autodata/protein_encoding/{}".format(seqfile), header=0, index_col=0)
    # add_dataframe["Entry"]=add_dataframe.index
    # add_dataframe.reset_index(drop=True,inplace=True)
    #drop all zero columns

    add_dataframe = (add_dataframe.loc[:, add_dataframe.sum() != 0.0])
    add_dataframe.drop_duplicates(subset="Entry", inplace=True)

    add_dataframe["Entry"] = (add_dataframe["Entry"].astype("string"))
    #convert into int
    add_dataframe.iloc[:,:-3]=add_dataframe.iloc[:,:-3].astype("int64")

    print(add_dataframe.dtypes)
    #add_dataframe.to_csv("../autodata/protein_encoding/{}".format(seqfile))

    print("merging fingerprint and sequences encoding---------")
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    #print(input_dataframe)
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    print("merged inputdataframe after drop na")
    print(input_dataframe)
    print("saving input data.......")
    #
    input_dataframe.to_csv(
        "../autodata/input_data/input{}fg_bond{}_{}_k_mer.csv".format(str(num_bit),
                                                            str(radius),
                                                            group))
    return input_dataframe
def sepreate_input(auto="",file="",numbit:int=2048,bond:int=2):
    """

    :param auto:
    :param file:
    :param numbit:
    :param bond:
    :return:
    """
    input_dataall=pd.read_csv(file,header=0,index_col=0)
    seprate_dataset = input_dataall.groupby(by=["methyl_type"])
    for group in seprate_dataset.groups:
        sub_df = seprate_dataset.get_group(group)
        group = sub_df["methyl_type"].unique()
        print(group)
        sub_df.reset_index(drop=True,inplace=True)

        sub_df.to_csv("../{}data/group/{}_{}_{}_with_bitinfo.csv".format(auto,group,str(numbit),str(bond)))
def perform_cluster_based_on_substrate(file_directory="../autodata/fingerprint_bit128_radius3_all_data_drop_atom.csv"):
    input_data=pd.read_csv("{}".format(file_directory),header=0,index_col=0)

    # Define clustering setup
    def ClusterFps(fps, cutoff=0.2):
        from rdkit import DataStructs
        from rdkit.ML.Cluster import Butina

        # first generate the distance matrix:
        dists = []
        nfps = len(fps)
        for i in range(1, nfps):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend([1 - x for x in sims])

        # now cluster the data:
        cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
        return cs
    fingerprintlist=[]
    fg_df = (input_data.iloc[:,:256])
    fgs=[]
    print(fg_df.shape)
    for i,index in enumerate(fg_df.index):
        line = fg_df.loc[[index]]
        print(line)
        fgs.append("".join([str(i) for i in line]))
        fgs[i]=DataStructs.cDataStructs.CreateFromBitString(fgs[i])
    print(fgs)
    clusters = ClusterFps(fgs, cutoff=0.4)
    print(clusters)
def cluster_with_initial_centroids(X):

    from sklearn.cluster import KMeans
    centroid_idx = [0, 2]  # let data point 0 and 2 be our centroids
    centroids = X[centroid_idx, :]
    print(centroids)  # [[1. 0. 0.]
    # [0. 0. 1.]]

    kmeans = KMeans(n_clusters=2, init=centroids,
                    max_iter=1)  # just run one k-Means iteration so that the centroids are not updated

    kmeans.fit(X)
    kmeans.labels_

def visilize_bit_info_of_created_fingerprint(smile,radius,numbit,atom_index):
    substrate_molecular=Chem.MolFromSmiles(smile)
    atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
        substrate_molecular,
        radius,
        atom_index
    )
    atom_map = {}
    submol = Chem.PathToSubmol(substrate_molecular, atom_environment,
                               atomMap=atom_map)
    Chem.SanitizeMol(submol,
                     sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    rdmolops.SanitizeFlags.SANITIZE_NONE
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect( submol,
            radius,
         numbit,bitInfo=bi
        )
    print(fp)
    print(list(fp.GetOnBits()))
    Draw.ShowMol(submol,(600,600))
    #226 is the No.1 feature importance
    im=Draw.DrawMorganBit(submol, (226-numbit), bi)
    im.show()
def main():
    today = date.today()
    # dd/mm/YY
    d1 = today.strftime("%d_%m_%Y")


    smile="[27*][26C@@H]1[25O][24C@H]([23CH2][22O][21P](=[34O])([35O-])[20O][19P](=[36O])([37O-])[18O][17P](=[38O])([39O-])[16O][15CH2][14C@H]2[13O][12C@@H]([11n]3[10cH][8n+]([9CH3])[7c]4c(=[1O])[2nH][3c]([4NH2])[5n][6c]43)[42C@H]([43OH])[40C@@H]2[41OH])[31C@@H]([32O][33*])[28C@H]1[29O]"
    mol1= Chem.MolFromSmiles(smile)
    Draw.ShowMol(mol1,(600,600))
    #visilize_bit_info_of_created_fingerprint(smile, 3,128,18)
    # mo_del = Model_class()
    # sequence_data = pd.read_csv("../autodata/protein_encoding/all_k_mer_encoding_sepreate_without_align.csv",header=0,index_col=0)
    # sequence_data['methyl_type']=sequence_data["group"]
    # print(sequence_data)
    # sequence_data.drop("group", inplace=True,axis=1)
    # print(sequence_data.columns)
    # mo_del.run_PCA(sequence_data,y_label=sequence_data['methyl_type'],file_name="all_k_mer_encoding_sepreate_without_align")

    """separate the fingerprint based on methylation group"""
    #sepreate_input("auto","../autodata/fingerprint_bit128_radius3_all_data_drop_atom.csv",128,3)

    # for group in ["N","O","S","C"]:
    #     use_k_merencoding_for_create_input(x="../autodata/group/['{}']_128_3_with_bitinfo.csv".format(group),num_bit=128,radius=3,seqfile="{}_k_mer_encoding_without_align_26_08.csv".format(group),group=group)
    #perform_cluster_based_on_substrate()
    #mo_del.hierarchical_clustering(sequence_data)
    # data_with_site = pd.read_csv("../data/mannual_data.csv", header=0,
    #                              index_col=0)
    # with open('../data/methyl_site_dictionary', 'rb') as file1:
    #     diction_atom = dill.load(file1)
    # data_with_site = data_with_site.fillna(0)
    # mo_del.save_fingerprints_to_dataframe(data_with_site,diction_atom,128,3,True,"{}_manual_drop_duplicate_drop_atom_withtype_bond{}".format(str(128),str(3)))
    #

    #seq=sequences()
    #seq.group_seq_based_on_methylated_type()
    #seq.group_fg_based_on_methylated_type("data/input_dataframe_withoutstructure_dropatoms2048_drop_duplicate_drop_atom_withtype_bond2.csv",2048,2)
    groups=["N_seed","O_seed","S_seed","C_seed"]
    # for group in groups:
    #     parse_data.read_msa_and_encoding("{}".format(group))
    #sepreate_input("auto","../autodata/input_dataframe_withoutstructure_dropatoms128_drop_duplicate_drop_atom_withtype_bond3.csv",128,3)
    #create protein encoding
    # for group in groups:
    #     parse_data.read_msa_and_encoding("{}".format(group))


    # groups1 = ["S","C","O","N"]
    # for group in groups1:
    #     print(group)
    #     build_different_input(auto="auto",x="../autodata/group/['{}']_128_3_with_bitinfo.csv".format(group),num_bit=128,radius=3,seqfile="{}_seed_onehot_encoding.csv".format(group),group=group)

    '''
    for file in groups:
        input_dataframe = pd.read_csv(
            "../autodata/input_data/input128fg_dpna_bond3_{}_onehot_encoding.csv.csv".format(file), header=0, index_col=0)

        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
            input_dataframe)

        # mo_del.three_D_pca(X_train, y_train, "{}_128_2".format(file))
        # mo_del.run_PCA(X_train, y_train, "{}_128_2".format(file))
        X_train = X_train.drop(columns=["methyl_type"])
        X_test = X_test.drop(columns=["methyl_type"])
        y_train = y_train.drop(columns=["methyl_type"])
        y_test = y_test.drop(columns=["methyl_type"])
        # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
        #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
        model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                                 "sepreate_align_input128fg_bi_type_bond3_rf{}_{}_spreate_seed".format(
                                     d1, file), i=0)
    '''
    #parse_data.read_msa_and_encoding(file_name="uniprot_and_manual_align")
    # mo_del.group_by_site()
    #sepreate_input(file="../autodata/input_data/input128fg_dpna_bond2_6_seed_onehot_encoding.csv.csv", numbit = 128, bond= 3)
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


