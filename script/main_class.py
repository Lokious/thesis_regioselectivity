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
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from IPython.display import SVG
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
            try:
                add_dataframe=pd.read_csv("../{}data/protein_encoding/{}".format(auto,seqfile),header=0,index_col=0)
            except:
                print("create sequence encoding from msa......")
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

        sub_df.to_csv("../{}data/group/{}_{}_{}_with_bitinfo_19_09.csv".format(auto,group,str(numbit),str(bond)))

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

def visilize_bit_info_of_created_fingerprint(sub_smile,radius,numbit,atom_index,product_smile,mol_id):

    substrate_molecular=Chem.MolFromSmiles(sub_smile)
    #product_molecular = Chem.MolFromSmiles(product_smile)
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
    print(list(fp.GetOnBits()))


    #226 is the No.1 feature importance
    try:
        im=Draw.DrawMorganBit(submol, (226-numbit), bi)
        #im.show()
        im.save("fig/{}_morganBit.png".format(str(mol_id)+"atom"+str(atom_index)))
        #Draw.ShowMol(substrate_molecular, (600, 600))
        Draw.MolToFile(substrate_molecular,"fig/{}_substarte.png".format(str(mol_id)+"atom"+str(atom_index)),(600,600))
        print(product_smile)
        # Draw.ShowMol(product_molecular, (600, 600))
        #Draw.ShowMol(submol, (600, 600))
        Draw.MolToFile(submol,
                       "fig/{}_atom_environment.png".format(str(mol_id)+"atom"+str(atom_index)), (600, 600))
    except:
        print("this bit is not in atom fingerprint")

def show_bitinfo_for_input_data(file="",smile_file=""):
    input_datafrme= pd.read_csv(file,index_col=0,header=0)
    smile_df=pd.read_csv(smile_file,index_col=0,header=0)
    entry_list= input_datafrme["Entry"].unique()
    for id in smile_df.index:
        if smile_df.loc[id,"Entry"] in entry_list:
            sub_smile=smile_df.loc[id,"main_sub"]
            reactant_site = smile_df.loc[id, "reactant_site"]
            mol=Chem.MolFromSmiles(sub_smile)
            mol_id = id
            for atom in mol.GetAtoms():
                atom_index= atom.GetIdx()
                visilize_bit_info_of_created_fingerprint(sub_smile, radius=3, numbit=128,
                                                     atom_index=atom_index,product_smile=reactant_site,mol_id=mol_id)
def check_if_atom_environment_has_specific_bit(substrate_molecule, radius=3,
                                                 numbit=128,
                                                 atom_index=0,bit:int=226):
    atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
        substrate_molecule,
        radius,
        atom_index
    )
    atom_map = {}
    submol = Chem.PathToSubmol(substrate_molecule, atom_environment,
                               atomMap=atom_map)
    Chem.SanitizeMol(submol,
                     sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    rdmolops.SanitizeFlags.SANITIZE_NONE
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(submol,
                                               radius,
                                               numbit, bitInfo=bi
                                               )
    #print(fp)
    bit_fingerprint_mol = np.zeros(
        (0,),
        dtype=int)  # (one dimention, 0 is number of rows)
    DataStructs.ConvertToNumpyArray(fp, bit_fingerprint_mol)
    fg=list(fp.GetOnBits())
    #print(fg)
    if (bit-128) in fg:
        return True
    else:
        return False
def highlight_atom_with_226_in_represent_mol(smile_file):
    smile_df = pd.read_csv(smile_file, index_col=0, header=0)
    smile_df.dropna(inplace=True)

    # unique_sub=list((smile_df["reactant_site"].value_counts()).index)
    unique_sub=[]
    for i in smile_df.index:
        substrate = smile_df.loc[i, "main_sub"]
        if substrate not in unique_sub:
            unique_sub.append(substrate)
            try:
                sites = smile_df.loc[i, "reactant_site"]
                print(sites)
                if "," in sites:
                    site_list = sites.split(",")
                else:
                    site_list = [sites]
                sites = []
                for site in site_list:
                    sites.append(site.split(":")[1])
                #print(sites)
                mol = Chem.MolFromSmiles(substrate)
                high_light_index = []
                high_light_colour={}
                site_index = []
                for atom in mol.GetAtoms():
                    atom_index = atom.GetIdx()
                    print("index{}".format(atom_index))
                    Boolean_226=check_if_atom_environment_has_specific_bit(mol,radius=3,numbit=128,
                                                                   atom_index=atom_index,
                                                                   bit= 226)
                    if Boolean_226:
                        high_light_index.append(atom_index)
                        high_light_colour[atom_index]=(1,0,0)
                    isotope=str(atom.GetIsotope())
                    print(isotope)
                    if isotope in sites:
                        print(atom_index)
                        site_index.append(atom_index)
                #site with 226 is red methyl site is green is both then yellow
                change_colour_index=[]
                change_colour_colour={}
                for methyl_site in site_index:
                    if methyl_site in high_light_index:
                        high_light_index.remove(methyl_site)
                        del high_light_colour[methyl_site]
                        change_colour_index.append(methyl_site)
                        change_colour_colour[methyl_site]=(1,1,0)
                    else:
                        change_colour_index.append(methyl_site)
                        change_colour_colour[methyl_site]=(0,1,0)
                high_light_index=high_light_index+change_colour_index
                high_light_colour.update(change_colour_colour)
                print(high_light_index)
                print(high_light_colour)
                drawer = Draw.MolDraw2DCairo(800, 800)
                drawer.DrawMolecule(mol, highlightAtoms=high_light_index,
                                    highlightAtomColors=high_light_colour)

                drawer.FinishDrawing()

                png = drawer.GetDrawingText()

                # save png to file
                with open('hight_light_226/{}.png'.format(i), 'wb') as png_file:
                    png_file.write(png)
                #
                # img1= Draw.MolToImage(mol,highlightAtoms=high_light_index,highlightAtomColors=high_light_colour,size=(600,600))
                # img1.save("hight_light_226/{}.jpg".format(i))
                print(substrate)
            except:
                print(i)
                #print(substrate)
                print("cannot get mol from smile")
def main():
    today = date.today()
    # dd/mm/YY
    d1 = today.strftime("%d_%m_%Y")

    # drawer = Draw.MolDraw2DCairo(800, 800)
    # print("here1")
    # mol=Chem.MolFromSmiles("C(C)C")
    # drawer.DrawMolecule(mol, highlightAtoms=[0],
    #                                   highlightAtomColors={0:(0,0,1)})
    # drawer.FinishDrawing()
    # svg = drawer.GetDrawingText()
    #
    #
    # drawer.DrawMolecule(mol, highlightAtoms=[0],
    #                                   highlightAtomColors={0:(0,0,1)})
    # drawer.FinishDrawing()
    # svg = drawer.GetDrawingText().replace('svg:', '')
    # SVG(svg)
    highlight_atom_with_226_in_represent_mol(smile_file="../autodata/seq_smiles_all.csv")
    #show_bitinfo_for_input_data(file="../autodata/input_data/active_site/PF08241PF01795_bit_score11_coverage0.7_ACS_bit128_3_remove_redundant.csv", smile_file="../autodata/seq_smiles_all.csv")
    #build_different_input(auto="auto", x="../autodata/group/['N']_128_3_with_bitinfo.csv",num_bit = 128, radius:int = 3, seqfile= "6_seed_onehot_encoding.csv", group = "N")
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
    #     build_different_input(auto="auto",x="../autodata/group/['{}']_128_3_with_bitinfo_19_09.csv".format(group),num_bit=128,radius=3,seqfile="{}_seed_onehot_encoding_sepreate.csv".format(group),group=group)

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
    #sepreate_input("auto","../autodata/fingerprint/fingerprint_bit128_radius3_all_data_drop_atom_19_09.csv",128,3)
    '''
    #for active site encoding
    #O methyltransferase
    X=pd.read_csv("../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data.csv",header=0,index_col=0)
    add_dataframe=pd.read_csv("../autodata/protein_encoding/6_seed_onehot_encoding.csv_128fg.csv",header=0,index_col=0)

    print(add_dataframe)
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    print(input_dataframe)
    input_dataframe = input_dataframe.dropna(axis=0,how="any")
    print(input_dataframe)
    input_dataframe.to_csv("../autodata/input_data/active_site/6_seed_onehot_encoding_MACCSkey_no_same_sub.csv")
    mo_del = Model_class()
    print(input_dataframe)
    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        input_dataframe,group_column="main_sub")

    X_train = X_train.drop(columns=["methyl_type", "molecular_id","atom_index"])
    # save x test for further analysis result
    y_test.to_csv(
        "../autodata/model/6_seed_onehot_encoding_MACCSkey_no_same_sub_y_test.csv")
    X_test.to_csv(
        "../autodata/model/6_seed_onehot_encoding_MACCSkey_no_same_sub_X_test.csv")
    X_test = X_test.drop(columns=["methyl_type", "molecular_id","atom_index"])
    # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
    model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                             "6_seed_onehot_encoding_MACCSkey_no_same_sub", i=0)

    '''
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


