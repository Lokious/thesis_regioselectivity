#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
datafile should be put in the directory data/ ,this script should be under regioselectivity_prediction
includes Rhea-ec_2_1_1.tsv
"""
import dill
import sklearn.metrics


from molecular_class import Molecule, Reaction,main_substrate
import parse_data
#for data structure
import pandas as pd
import numpy as np

import random
from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
#import pikachu
from rdkit.Chem import Draw, AllChem, rdmolops
import glob
import joblib
import pathlib
#for display and draw picture and graph
from IPython.display import SVG, display, Image
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

from PIL import Image
import copy
#import for model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import confusion_matrix, mean_squared_error, mean_absolute_percentage_error, r2_score, mean_absolute_error, accuracy_score, ConfusionMatrixDisplay, multilabel_confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import GroupShuffleSplit
from sklearn.decomposition import PCA
from sklearn.metrics import make_scorer
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import precision_score, accuracy_score, plot_roc_curve, RocCurveDisplay,ConfusionMatrixDisplay, roc_auc_score, roc_curve
from sklearn.svm import SVC
from sklearn.svm import SVC
import matplotlib.patches as mpatches
from math import sqrt
import unittest

class Model_class():

    #def __init__(self):
        #self.check_file_exist()

    def check_file_exist(self,file1='../data/seq_smiles_all_MANUAL.csv',file2='../data/diction_atom_all'):
        """
        This function is used to check if the file exit, otherwise create the files

        :param file1: path of csv file, which includes substrates and sequneces information
        :param file2: object file, which is a dictionary includes index from
        file1 as key and atom object(which is methylated in substrate molecular) as value
        :return: None
        """


        try:
            assert pathlib.Path(file1).exists()==True
            assert pathlib.Path(file2).exists()==True

        except:
            print("missing file{},{}".format(file1,file2))
            print("preparing the file......")
            rh_file = "../autodata/rawdata/rhea2uniprot_sprot.tsv"
            rheauniprot_dataframe = parse_data.readrhlist(rh_file)
            #read id and sequences in dataframe
            seq_file = "../autodata/rawdata/id_tosequence.xlsx" # run with pycharm
            id_seq_dataframe = parse_data.read_sequence(seq_file)
            seq_smiles = self.merge_uniprot_id_smile(rheauniprot_dataframe,id_seq_dataframe)
            #data_frame = self.keep_longest_smile(seq_smiles)
            #data_frame=self.keep_methyled_substrate(seq_smiles)
            seq_smiles.to_csv("../autodata/seq_smiles_all.csv")
            self.return_reactions(seq_smiles)



    def merge_uniprot_id_smile(self,rheauniprot_dataframe,seq_df):
        """
        combine sequence and substrate smile in onefile and save to csv file
        :param
            rheauniprot_dataframe:pd dataframe includes uniprot entry and rhea id
            seq_df: dataframe which includes sequences and uniprot entry
        :return:
        """
        df1 = parse_data.get_substrate_chebi("../data/rawdata/Rhea-ec_2_1_1.tsv")
        df1 = parse_data.get_smile(df1)
        df1.index = df1.index.map(lambda x: x.split(":")[1])
        df1["RHEA_ID"] = df1.index
        df1.set_index("RHEA_ID", inplace=True)
        #list of uniprot entry
        uniprot_entry=pd.read_csv("../autodata/rawdata/uniprot-ec2.1.1.tsv",header=0,sep="\t")["Entry"].tolist()
        #print(uniprot_entry)
        # uniprot_entry = parse_data.remove_duplicated_id(
        #     r"E:\Download\regioselectivity_prediction\data\hmm_out")
        fulldata = parse_data.rheaid_to_uniprot(uniprot_entry,
                                                rheauniprot_dataframe)
        fulldata.set_index("RHEA_ID", inplace=True)
        fulldata = fulldata.reset_index()
        df1 = df1.reset_index()
        complete_data = pd.merge(fulldata, df1)
        seq_df = seq_df.reset_index()
        complete_data = pd.merge(seq_df,complete_data)
        complete_data = parse_data.merge_reaction(complete_data)

        return complete_data

    def keep_longest_smile(self,dataframe_before):

        raise RuntimeError(
            "function `keep_longest_smile()` is deprecated")
        dataframe_before["main_sub"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        dataframe_before["main_pro"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        for index in dataframe_before.index:
            main_sub = max((dataframe_before.loc[index,"sub_smiles"]).split("."), key=len)
            dataframe_before.loc[index, "main_sub"] = main_sub
            main_pro = max((dataframe_before.loc[index,"pro_smiles"]).split("."), key=len)
            dataframe_before.loc[index, "main_pro"] = main_pro
        return dataframe_before

    def keep_methyled_substrate(self,dataframe_before):
        """
        exclude H+,keep the substrate which is methylated not the methyl donor

        :param dataframe_before: substrate includes mainsubstrates
        :return:
        """
        raise RuntimeError("function `keep_methyled_substrate()` is deprecated")
        dataframe_before["main_sub"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        dataframe_before["main_pro"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        for index in dataframe_before.index:
            subs = (dataframe_before.loc[index,"sub_smiles"]).split(".")
            pros = (dataframe_before.loc[index,"pro_smiles"]).split(".")
            for sub in subs:
                mol = Chem.MolFromSmiles(sub)
                #remove very small molecular
                if len(mol.GetAtoms())<2:
                    subs.remove(sub)
            for pro in pros:
                mol = Chem.MolFromSmiles(pro)
                if len(mol.GetAtoms())<2:
                    pros.remove(pro)
            if len(subs) == len(pros):
                #return the smile in list
                smiles = main_substrate(subs, pros)
                if smiles:
                    main_sub, main_pro = smiles[0],smiles[1]
                    dataframe_before.loc[index, "main_sub"] = main_sub
                    dataframe_before.loc[index, "main_pro"] = main_pro

        return copy.deepcopy(dataframe_before)

    def return_reactions(self,dataframe_rr,test=False):
        """
        This function is to save methylated site in file

        :param dataframe:
        :return:None
        """
        atom_object_dictionary = {}
        dataframe_rr["reactant_site"] = pd.DataFrame(
            len(dataframe_rr.index) * [0]).astype('object')
        # dataframe_rr["mainsub_mol"] = pd.DataFrame(
        #     len(dataframe_rr.index) * [0]).astype('object')
        for index in dataframe_rr.index:

            subs = dataframe_rr.loc[index, "sub_smiles"]
            pros = dataframe_rr.loc[index, "pro_smiles"]
            if (subs != pros) and len(subs) !=0 and len(pros) !=0:
                reaction1 = Reaction()
                #print(dataframe.loc[index,"RHEA_ID"])
                #rxn_file_name = "../data/rxn_picture/{}".format(dataframe_rr.loc[index,"Entry"])
                #r1 = reaction1.get_reaction_sites(rxn_object=rxn,file_name=rxn_file_name)
                print("index: {}".format(index))
                pro_mol,remove_methyl_smile,list_methylsite,check = reaction1.get_reaction_sites(pros,subs)
                Draw.MolToFile(pro_mol,
                               "../pro_fig/{}_pro.png".format(index))
                rxn = AllChem.ReactionFromSmarts((remove_methyl_smile+">>"+Chem.MolToSmiles(pro_mol)))
                d2d=Draw.MolDraw2DCairo(800, 800)
                d2d.DrawReaction(rxn)
                png = d2d.GetDrawingText()
                open('../autodata/rxn_picture/reaction{}.png'.format(index), 'wb+').write(png)
                #Draw.MolToFile(mainsub_mol,"../data/substrate_mol_picture/{}.png".format(index),size=(600,600))
                #save atom index(methylation site) from substrate in dataframe
                site = ",".join(list_methylsite)
                #print(site)
                #dataframe_rr.loc[index,"mainsub_mol"] = mainsub_mol
                dataframe_rr.loc[index, "main_sub"] =remove_methyl_smile
                dataframe_rr.loc[index, "check"] = check
            else:
                site = ""
            if site != "":
                dataframe_rr.loc[index, "reactant_site"] = site
                import re
                methyl_type = set((re.sub(r'[^a-zA-Z]', ' ', site)).split())
                print(methyl_type)
                dataframe_rr.loc[index, "methyl_type"]="_".join(methyl_type)
            else:
                #because we use the largest molecular, but for some substrates, methyl donor is larger
                #we first leave those
                dataframe_rr.loc[index,"reactant_site"] = "NA"
            #save list of methyl site atom objects and index to a dictionary
            atom_object_dictionary[index] = site
        else:
            if test == False:
                with open("../autodata/seq_smiles_all", "wb") as dill_file:
                                    dill.dump(dataframe_rr, dill_file)
                with open("../autodata/diction_atom_all", "wb") as dill_file:
                    dill.dump(atom_object_dictionary, dill_file)
                dataframe_rr.to_csv("../autodata/seq_smiles_all.csv")
            else:
                print("won't save file while running test")
        return dataframe_rr

    def save_fingerprints_to_dataframe(self,sauce_data,atom_object_dictionary,num_bits: int = 2048,radius: int = 3,drop_atoms=False,file_name=""):
        """
        this function is to build inputdata with fingerprints and labels

        :param sauce_data: pd.dataframe with substrate and metylation site information
        :param atom_object_dictionary: key: index of sauce_data value: mathylation site
        :param num_bits: num_bits for morgen fingerprint
        :param radius:
        :return:
        """

        self_defined_mol_object = Molecule()
        input_dataframe = pd.DataFrame()
        current_index = 0

        for index in sauce_data.index:

            print(index)
            try:
                sub_mol = Chem.MolFromSmiles(sauce_data.loc[index,"main_sub"])
                # print(sub_mol)
                #Draw.ShowMol(sub_mol, size=(600, 600))
                sub_rest_mol, no_use_variable = self_defined_mol_object.mol_with_atom_index(copy.deepcopy(sub_mol))
                fingerprint_mol = self_defined_mol_object.create_fingerprint_mol(
                    sub_rest_mol, num_bits=num_bits, radius=radius)
                for atom in sub_mol.GetAtoms():
                    #set label
                    #isotope is what we saved before in the list of methylaion site
                    sy_index =(atom.GetSymbol() + ":" + str(atom.GetIsotope()))
                    # print(sy_index)
                    # print(atom_object_dictionary[index])
                    if sy_index in atom_object_dictionary[index]:
                        print(sy_index)
                        label = 1
                    else:
                        label = 0
                    #atom_index_sub = atom.GetAtomMapNum()
                    newrow = {}
                    #if drop_atoms is True, filter out some Carbon whose degree is already 4
                    if drop_atoms and (label != 1):
                        #print(atom.GetTotalDegree(),atom.GetSymbol())
                        if (atom.GetTotalDegree() == 4) and (atom.GetSymbol()=="C"):
                            print("drop C")
                            continue
                        else:
                            # resrt atom index and then build fingerprint
                            fingerprint_atom = self_defined_mol_object.create_fingerprint_atom(
                                sub_rest_mol, atom_object=atom, num_bits=num_bits,
                                radius=radius)
                            # add fingerprint atom index ebedding sequences and label to dataframe
                            for i, item in enumerate(fingerprint_mol):
                                newrow[i] = item

                            for j, item in enumerate(fingerprint_atom):
                                newrow[j + i+1] = item

                            # newrow['atom_index'] = atom_index_sub
                            add_dataframe = pd.DataFrame(newrow,
                                                         index=[current_index])

                            input_dataframe = pd.concat(
                                [input_dataframe, add_dataframe], axis=
                                0)
                            input_dataframe.loc[
                                current_index, "molecular_id"] = "m" + str(index)
                            input_dataframe.loc[current_index, "label"] = label
                            input_dataframe.loc[current_index, "Entry"] = \
                            sauce_data.loc[index, "Entry"]
                            input_dataframe.loc[current_index, "methyl_type"] = sauce_data.loc[index, "methyl_type"]
                            input_dataframe.loc[
                                current_index, 'atom_index'] = sy_index
                            current_index += 1

                    else:
                        #resrt atom index and then build fingerprint
                        fingerprint_atom = self_defined_mol_object.create_fingerprint_atom(
                            sub_rest_mol, atom_object=atom, num_bits=num_bits, radius=radius)
                        #add fingerprint atom index ebedding sequences and label to dataframe
                        for i,item in enumerate(fingerprint_mol):
                            newrow[i] = item

                        for j,item in enumerate(fingerprint_atom):
                            newrow[j+i+1] = item

                        #newrow['atom_index'] = atom_index_sub
                        add_dataframe = pd.DataFrame(newrow,index=[current_index])

                        input_dataframe=pd.concat([input_dataframe,add_dataframe],axis=
                                                  0)
                        input_dataframe.loc[current_index, "molecular_id"] = "m"+str(index)
                        input_dataframe.loc[current_index,"label"] = label
                        input_dataframe.loc[current_index,"Entry"] = sauce_data.loc[index,"Entry"]
                        input_dataframe.loc[current_index, "methyl_type"] = sauce_data.loc[index, "methyl_type"]
                        input_dataframe.loc[
                            current_index, 'atom_index'] = sy_index
                        current_index += 1
                        print(current_index)
            except:
                    print("somethingwrong with this index{}".format(index))
                    continue
        if drop_atoms:
            input_dataframe.to_csv("../autodata/fingerprint/fingerprint_bit{}_radius{}_{}.csv".format(num_bits,radius,file_name))
            with open("../autodata/fingerprint_bit{}_radius{}_{}".format(num_bits,radius,file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        else:
            input_dataframe.to_csv("../autodata/fingerprint/fingerprint_bit{}_radius{}_{}.csv".format(num_bits,radius,file_name))
            with open("../autodata/fingerprint_bit{}_radius{}_{}".format(num_bits,radius,file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        return input_dataframe

    def save_MACCSkeys_fingerprints_to_dataframe(self,sauce_data,atom_object_dictionary,radius: int = 3,drop_atoms=False,file_name=""):
        """
        this function is to build inputdata with fingerprints and labels

        :param sauce_data: pd.dataframe with substrate and metylation site information
        :param atom_object_dictionary: key: index of sauce_data value: mathylation site
        :param num_bits: num_bits for morgen fingerprint
        :param radius:
        :return:
        """

        self_defined_mol_object = Molecule()
        input_dataframe = pd.DataFrame()
        current_index = 0

        for index in sauce_data.index:

            print(index)
            try:
                sub_mol = Chem.MolFromSmiles(sauce_data.loc[index,"main_sub"])
                # print(sub_mol)
                #Draw.ShowMol(sub_mol, size=(600, 600))
                sub_rest_mol, no_use_variable = self_defined_mol_object.mol_with_atom_index(copy.deepcopy(sub_mol))
                fingerprint_mol = self_defined_mol_object.create_MACCSkey_mol(
                    sub_rest_mol)
                for atom in sub_mol.GetAtoms():
                    #set label
                    #isotope is what we saved before in the list of methylaion site
                    sy_index = (atom.GetSymbol() + ":" + str(atom.GetIsotope()))
                    # print(sy_index)
                    # print(atom_object_dictionary[index])
                    if sy_index in atom_object_dictionary[index]:
                        print(sy_index)
                        label = 1
                    else:
                        label = 0
                    #atom_index_sub = atom.GetAtomMapNum()
                    newrow = {}
                    #if drop_atoms is True, filter out some Carbon whose degree is already 4
                    if drop_atoms and (label != 1):
                        #print(atom.GetTotalDegree(),atom.GetSymbol())
                        if (atom.GetTotalDegree() == 4) and (atom.GetSymbol()=="C"):
                            print("drop C")
                            continue
                        else:
                            # reset atom index and then build fingerprint
                            fingerprint_atom = self_defined_mol_object.create_MACCSkey_atom(
                                sub_rest_mol, atom_object=atom,
                                radius=radius)
                            # add fingerprint atom index ebedding sequences and label to dataframe
                            for i, item in enumerate(fingerprint_mol):
                                newrow[i] = item

                            for j, item in enumerate(fingerprint_atom):
                                newrow[j + i+1] = item

                            # newrow['atom_index'] = atom_index_sub
                            add_dataframe = pd.DataFrame(newrow,
                                                         index=[current_index])

                            input_dataframe = pd.concat(
                                [input_dataframe, add_dataframe], axis=
                                0)
                            input_dataframe.loc[
                                current_index, "molecular_id"] = "m" + str(index)
                            input_dataframe.loc[current_index, "label"] = label
                            input_dataframe.loc[current_index, "Entry"] = \
                            sauce_data.loc[index, "Entry"]
                            input_dataframe.loc[current_index, "methyl_type"] = sauce_data.loc[index, "methyl_type"]
                            input_dataframe.loc[
                                current_index, 'atom_index'] = sy_index
                            current_index += 1

                    else:
                        #reset atom index and then build fingerprint
                        fingerprint_atom = self_defined_mol_object.create_MACCSkey_atom(
                            sub_rest_mol, atom_object=atom,radius=radius)
                        #add fingerprint atom index ebedding sequences and label to dataframe
                        for i,item in enumerate(fingerprint_mol):
                            newrow[i] = item

                        for j,item in enumerate(fingerprint_atom):
                            newrow[j+i+1] = item

                        #newrow['atom_index'] = atom_index_sub
                        add_dataframe = pd.DataFrame(newrow,index=[current_index])

                        input_dataframe=pd.concat([input_dataframe,add_dataframe],axis=
                                                  0)
                        input_dataframe.loc[current_index, "molecular_id"] = "m"+str(index)
                        input_dataframe.loc[current_index,"label"] = label
                        input_dataframe.loc[current_index,"Entry"] = sauce_data.loc[index,"Entry"]
                        input_dataframe.loc[current_index, "methyl_type"] = sauce_data.loc[index, "methyl_type"]
                        input_dataframe.loc[
                            current_index, 'atom_index'] = sy_index
                        current_index += 1
                        print(current_index)
            except:
                    print("somethingwrong with this index{}".format(index))
                    continue
        if drop_atoms:
            input_dataframe.to_csv("../autodata/fingerprint/MACCS_fingerprint_bit{}_radius{}_{}_2_11.csv".format(167,radius,file_name))
            with open("../autodata/MACCS_fingerprint_bit{}_radius{}_{}_2_11".format(167,radius,file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        else:
            input_dataframe.to_csv("../autodata/fingerprint/MACCS_fingerprint_bit{}_radius{}_{}_2_11.csv".format(167,radius,file_name))
            with open("../autodata/MACCS_fingerprint_bit{}_radius{}_{}_2_11".format(167,radius,file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        return input_dataframe

    def fingerprint_df_preapare(self,substrate_mol_df,num_bits: int = 2048,radius: int = 3):
        """
        This function is use for create fingerprint for given dataframe with substrate smile and reactant site
        not use in prepare training and test data
        :param substrate_mol:
        :param num_bits:
        :param radius:
        :return:
        """
        current_index = 0
        mol_id:int = 0
        fingerprint_dataframe = pd.DataFrame()
        self_defined_mol_object = Molecule()
        print(substrate_mol_df)
        for index in substrate_mol_df.index:
            substrate_mol = substrate_mol_df.loc[index,'substrate_mol']
            try:

                #Draw.ShowMol(substrate_mol, size=(600, 600))
                fingerprint_mol = self_defined_mol_object.create_fingerprint_mol(
                    substrate_mol, num_bits=num_bits, radius=radius)
                for atom in substrate_mol.GetAtoms():
                    try:
                        sites=substrate_mol_df.loc[index,'reactant_site'].split(";")
                        methyl_sites=[(x.split(":")[0][0] + ":" +x.split(":")[1][0:]) for x in sites]
                        print(methyl_sites)
                    except:
                        symbol_index=substrate_mol_df.loc[
                            index, 'reactant_site']
                        methyl_sites = [symbol_index.split(":")[0][0] + ":" +symbol_index.split(":")[1][0:]]
                    sy_index = (atom.GetSymbol() + ":" + str(atom.GetIsotope()))
                    #print(sy_index)
                    if sy_index in methyl_sites:
                        print(sy_index)
                        label = 1
                    else:
                        label = 0
                    fingerprint_atom = self_defined_mol_object.create_fingerprint_atom(
                        substrate_mol, atom_object=atom, num_bits=num_bits, radius=radius)
                    atom_index_sub = atom.GetIdx()
                    newrow = {}
                    # add fingerprint atom index ebedding sequences and label to dataframe
                    for i, item in enumerate(fingerprint_mol):
                        newrow[i] = item
                    for j, item in enumerate(fingerprint_atom):
                        newrow[j + i+1] = item

                    add_dataframe = pd.DataFrame(newrow, index=[current_index])
                    fingerprint_dataframe = pd.concat([fingerprint_dataframe, add_dataframe], axis=
                    0)
                    fingerprint_dataframe.loc[current_index, "molecular_id"] = "m"+str(mol_id)
                    fingerprint_dataframe.loc[current_index, 'atom_index'] = atom_index_sub
                    fingerprint_dataframe.loc[
                        current_index, 'label'] = label
                    current_index += 1
            except:
                print("skip in mol{}".format(mol_id))
            mol_id +=1
        print(fingerprint_dataframe)
        return fingerprint_dataframe

    def MACCSkey_fingerprint_df_preapare(self,substrate_mol_df,radius: int = 3):
        """
        This function is use for create fingerprint for given dataframe with substrate smile and reactant site
        not use in prepare training and test data
        :param substrate_mol:
        :param num_bits:
        :param radius:
        :return:
        """
        current_index = 0
        mol_id:int = 0
        fingerprint_dataframe = pd.DataFrame()
        self_defined_mol_object = Molecule()
        print(substrate_mol_df)
        for index in substrate_mol_df.index:
            substrate_mol = substrate_mol_df.loc[index,'substrate_mol']
            try:

                #Draw.ShowMol(substrate_mol, size=(600, 600))
                fingerprint_mol = self_defined_mol_object.create_MACCSkey_mol(
                    substrate_mol)
                for atom in substrate_mol.GetAtoms():
                    try:
                        sites=substrate_mol_df.loc[index,'reactant_site'].split(";")
                        methyl_sites=[(x.split(":")[0][0] + ":" +x.split(":")[1][0:]) for x in sites]
                        print(methyl_sites)
                    except:
                        symbol_index=substrate_mol_df.loc[
                            index, 'reactant_site']
                        methyl_sites = [symbol_index.split(":")[0][0] + ":" +symbol_index.split(":")[1][0:]]
                    sy_index = (atom.GetSymbol() + ":" + str(atom.GetIsotope()))
                    #print(sy_index)
                    if sy_index in methyl_sites:
                        print(sy_index)
                        label = 1
                    else:
                        label = 0
                    fingerprint_atom = self_defined_mol_object.create_MACCSkey_atom(
                        substrate_mol, atom_object=atom, radius=radius)
                    atom_index_sub = atom.GetIdx()
                    newrow = {}
                    # add fingerprint atom index ebedding sequences and label to dataframe
                    for i, item in enumerate(fingerprint_mol):
                        newrow[i] = item
                    for j, item in enumerate(fingerprint_atom):
                        newrow[j + i+1] = item

                    add_dataframe = pd.DataFrame(newrow, index=[current_index])
                    fingerprint_dataframe = pd.concat([fingerprint_dataframe, add_dataframe], axis=
                    0)
                    fingerprint_dataframe.loc[current_index, "molecular_id"] = "m"+str(mol_id)
                    fingerprint_dataframe.loc[current_index, 'atom_index'] = atom_index_sub
                    fingerprint_dataframe.loc[
                        current_index, 'label'] = label
                    current_index += 1
            except:
                print("skip in mol{}".format(mol_id))
            mol_id +=1
        print(fingerprint_dataframe)
        return fingerprint_dataframe


    def prepare_train_teat_data(self,df,test:float=0.2,group_column:str='molecular_id',i:int=0):
        """
        simply sepreate train and test data
        :param data:
        :return:
        """

        if group_column=="molecular_id":
            splitter = GroupShuffleSplit(test_size=test, n_splits=1, random_state=i)
            split = splitter.split(df, groups=df[group_column])
            train_inds, test_inds = next(split)
            train = df.iloc[train_inds]
            test = df.iloc[test_inds]

            #X_train = train[list(train.columns)[:-2]]
            X_train = (copy.deepcopy(train)).drop(columns=["Entry","label"])
            Y_train = train["label"]
            #X_test = test[list(test.columns)[:-2]]
            X_test = (copy.deepcopy(test)).drop(columns=["Entry", "label"])
            Y_test = test["label"]
            #print(X_train, X_test, Y_train, Y_test)
            return X_train, X_test, Y_train, Y_test

        elif group_column=="main_sub":
            seq_smile_df = pd.read_csv("../autodata/seq_smiles_all.csv",header=0,index_col=0)
            seq_smile_df = pd.DataFrame(seq_smile_df,columns=["main_sub"])
            # print(seq_smile_df)
            # print(df)

            df["main_sub"]=pd.DataFrame(
                len(df.index) * [""]).astype('string')
            for i in df.index:
                index=int("".join(list(df.loc[i,"molecular_id"])[1:]))
                df.loc[i,"main_sub"]=seq_smile_df.loc[index,"main_sub"]
            #print(df)

            splitter = GroupShuffleSplit(test_size=test, n_splits=1,
                                         random_state=i)
            split = splitter.split(df, groups=df[group_column])
            train_inds, test_inds = next(split)
            train = df.iloc[train_inds]
            test = df.iloc[test_inds]
            print("number of substrates in train".format(len(train["main_sub"].unique())))
            print(len(test["main_sub"].unique()))
            X_train = (copy.deepcopy(train)).drop(columns=["Entry", "label","main_sub"])
            Y_train = train["label"]
            X_test = (copy.deepcopy(test)).drop(columns=["Entry", "label","main_sub"])
            Y_test = test["label"]
            return X_train, X_test, Y_train, Y_test
    def _get_colors(self,num_colors):
        import colorsys
        colors = []
        random.seed(0)
        for i in np.arange(0., 360., 360. / num_colors):
            hue = i / 360.
            lightness = (50 + np.random.rand() * 10) / 100.
            saturation = (90 + np.random.rand() * 10) / 100.
            colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        return colors

    def run_PCA(self,datasets,y_label,file_name=""):
        """

        :param datasets:
        :return:
        """
        #show pca result based on methylation type
        data_with_site = copy.deepcopy(datasets).drop(
            columns=["methyl_type"])
        def label_with_methylation_type():

            map_dictionary = {}

            colours = self._get_colors(len(datasets["methyl_type"].unique()))
            for type in datasets["methyl_type"].unique():
                map_dictionary[type] = colours.pop(0)
            colour_label = datasets["methyl_type"].map(map_dictionary)
            print(colour_label)
            return colour_label,map_dictionary

        def label_with_label():

            map_dictionary = {}
            label_df = copy.deepcopy(y_label)
            colours = self._get_colors(len(label_df.unique()))
            for type in label_df.unique():
                map_dictionary[type] = colours.pop(0)
            colour_label = label_df.map(map_dictionary)

            return colour_label,map_dictionary
        colour_label, map_dictionary=label_with_methylation_type()
        V = []
        PC = []
        #the most is 522 for pca

        for i in range(len(data_with_site.columns)):
            PC.append("PC" + str(i + 1))
            V.append("V" + str(i + 1))
            if i == min(521,(len(datasets.index)-1)):
                break

        pca_fit = PCA(n_components=len(PC),random_state=42).fit(data_with_site)

        pca_loadings = pd.DataFrame(pca_fit.components_.T,
                                    index=data_with_site.columns, columns=V)
        #show the lodings
        V1_PCA = pd.DataFrame(data=pca_loadings.abs(),
                              index=pca_loadings.index,
                              columns=["V1","V2"])

        V1_PCA['V_sum'] = V1_PCA.apply(lambda x: (x.abs()).sum(), axis=1)
        V1_PCA = V1_PCA.sort_values(by=["V_sum"], ascending=False)
        sns.barplot(data=V1_PCA.iloc[:30], x="V_sum", y=V1_PCA.index[:30]).set(
            title='Sum of v1 and v2 loading value of PCA')
        plt.tight_layout()
        plt.savefig("../Sum_of_v1_v2_loading_value_of_PCA_{}".format(file_name))

        plt.close()
        pca_df = pd.DataFrame(pca_fit.fit_transform(data_with_site),index=data_with_site.index, columns=PC)

        fig, ax = plt.subplots()
        plt.scatter(pca_df.PC1, pca_df.PC2, c=colour_label, s=3)

        #show explained variance
        handle_list = []
        for key in map_dictionary.keys():
            handle_list.append(mpatches.Patch(color=map_dictionary[key], label=key))
        ax.legend(
            handles=handle_list,loc="upper right",title="Type")

        #plt.scatter(pca_df.PC1, pca_df.PC2,  s=5,c=colour_label)
        print(pca_fit.explained_variance_ratio_)
        plt.xlabel("pc1({:.2f}%)".format(round(pca_fit.explained_variance_ratio_[0]*100),2))
        plt.ylabel("pc2({:.2f}%)".format(round(pca_fit.explained_variance_ratio_[1]*100),2))
        plt.title(" PCA coloured by methylation type")
        plt.savefig(
            "../pca_for encoding sequences and fingerprint_PC1 PC2{}".format(file_name))
        plt.clf()
        plt.close()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 pca_fit.explained_variance_ratio_, '-ro')
        plt.ylabel('Proportion of Variance Explained_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            '../Proportion of Variance Explained_{}'.format(file_name))

        plt.clf()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 np.cumsum(pca_fit.explained_variance_ratio_), '-o')
        plt.ylabel('Cumulative Proportion of Variance Explained_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            '../Cumulative Proportion of Variance Explained_{}'.format(file_name))

        plt.clf()
        ##pac with label
        colour_label, map_dictionary = label_with_label()
        V = []
        PC = []
        # the most is 522 for pca
        for i in range(len(data_with_site.columns)):
            PC.append("PC" + str(i + 1))
            V.append("V" + str(i + 1))
            if i == min(521,(len(datasets.index)-1)):
                break

        pca_fit = PCA(n_components=len(PC), random_state=42).fit(
            data_with_site)

        pca_df = pd.DataFrame(pca_fit.fit_transform(data_with_site),
                              index=data_with_site.index, columns=PC)
        fig, ax = plt.subplots()
        plt.scatter(pca_df.PC1, pca_df.PC2, c=colour_label, s=5)
        handle_list = []
        for key in map_dictionary.keys():
            handle_list.append(
                mpatches.Patch(color=map_dictionary[key], label=key))
        ax.legend(
            handles=handle_list)

        # plt.scatter(pca_df.PC1, pca_df.PC2,  s=5,c=colour_label)
        plt.xlabel("pc1({:.2f}%)".format(round(pca_fit.explained_variance_ratio_[0]*100),2))
        plt.ylabel("pc2({:.2f}%)".format(round(pca_fit.explained_variance_ratio_[1]*100),2))
        plt.title("First two component of PCA coloured by label type")
        plt.savefig(
            "../pca_for encoding sequences and fingerprint for label_{}".format(file_name))
        #plt.show()
        plt.close()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 pca_fit.explained_variance_ratio_, '-ro')
        plt.ylabel('Proportion of Variance Explained for label_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            '../Proportion of Variance Explained_{}'.format(file_name))
        #plt.show()
        plt.clf()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 np.cumsum(pca_fit.explained_variance_ratio_), '-o')
        plt.ylabel(
            '../Cumulative Proportion of Variance Explained foe label_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            '../Cumulative Proportion of Variance Explained foe label_{}'.format(file_name))
        #plt.show()
        plt.close()

        return pca_df

    def three_D_pca(self,datasets,y_label,file_name=""):

        import plotly.express as px
        from sklearn.decomposition import PCA


        X = copy.deepcopy(datasets).drop(
            columns=["methyl_type"])
        pca = PCA(n_components=3)
        components = pca.fit_transform(X)


        total_var = pca.explained_variance_ratio_.sum() * 100

        fig = px.scatter_3d(
            components, x=0, y=1, z=2, color=datasets["methyl_type"],
            title=f'Total Explained Variance: {total_var:.2f}%',
            labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
        )
        #fig.show()
        fig = px.scatter_3d(
            components, x=0, y=1, z=2, color=y_label,
            title=f'Total Explained Variance: {total_var:.2f}%',
            labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
        )
        #fig.show()



    def RF_model(self,X_train, X_test, y_train, y_test,file_name="",i:int=0):
        """

        :param X_train:
        :param X_test:
        :param y_train:
        :param y_test:
        :return: randomforest model
        """

        hyperparameters = {'n_estimators': [500,1000,1500,2000],
                           'max_features': [0.3,0.5,0.7],
                           'max_depth' : [10,15,20,30]
                           }

        #just for test
        # hyperparameters = {'n_estimators': [100,200],
        #                    'max_features': [0.3,0.5],
        #                    'max_depth' : [6,8]
        #                    }

        #n_jobs number of cores used for it
        #scoring is the Strategy to evaluate the performance of the cross-validated model on the test set
        rf_cv = GridSearchCV(RandomForestClassifier(random_state=0,class_weight="balanced",),
                             hyperparameters, scoring='roc_auc',
                             cv=5,
                             verbose=3,
                             n_jobs=12)

        ####use MCC as scoring#####
        # rf_cv = GridSearchCV(RandomForestClassifier(random_state=0,class_weight="balanced"),
        #                      hyperparameters, scoring=make_scorer(matthews_corrcoef),
        #                      cv=3,
        #                      verbose=3,
        #                      n_jobs=14)
        rf_cv.fit(X_train, y_train)
        print(rf_cv.best_params_)
        # roc for train data
        threshold_train = self.show_roc(rf_cv.best_estimator_, X_train, y_train, (file_name+"train"))
        #roc for test data
        threshold_test = self.show_roc(rf_cv.best_estimator_, X_test, y_test, (file_name + "test"))

        #confusion matrix for train
        self.cm_threshold(threshold_train, X_train, y_train, rf_cv.best_estimator_, (file_name+"train"))
        # self.cm_threshold(0.5, X_train, y_train, rf_cv.best_estimator_,
        #              (file_name + "train"))
        self.cm_threshold(threshold_test, X_test, y_test, rf_cv.best_estimator_, (file_name+"test"))
        cm_matrix=self.cm_threshold(threshold_train, X_test, y_test, rf_cv.best_estimator_,
                     (file_name + "test (same as train threshold)"))


        print(rf_cv.cv_results_)
        sns.lineplot(y=rf_cv.cv_results_["mean_test_score"],
                     x=rf_cv.cv_results_['param_max_depth'].data,
                     hue=rf_cv.cv_results_['param_n_estimators'])
        plt.xlabel("max_depth")
        plt.ylabel("roc_auc (mean 3-fold CV)")
        plt.title(
            "roc_auc score with different max_depth and features for RF model")
        plt.savefig("../autodata/separate_seed_result/Accuracy with different estimators and features for RF model_{}.svg".format(file_name))
        plt.close()

        #plt.show()

        #feature importance
        fi = pd.DataFrame(data=rf_cv.best_estimator_.feature_importances_, index=X_train.columns,
                          columns=['Importance']) \
            .sort_values(by=['Importance'], ascending=False)

        #And visualize
        ### size and font size should be set before the sns plot###
        plt.figure(figsize=(50,30))
        #plt.rcParams["figure.figsize"] = (40, 20)
        # sns_default=sns.axes_style()
        # print(sns_default)

        #sns.set(rc={'figure.figsize': (40,20)})
        sns.set(font_scale=3)
        ax1=sns.barplot(data=fi.head(30), x="Importance", y=(fi.head(30)).index)
        ax1.set_title(
            "feature importance for RF model")
        plt.tight_layout()
        #fig = ax.get_figure()

        plt.savefig('../autodata/separate_seed_result/feature_importance{}.svg'.format(file_name),bbox_inches = 'tight')
        plt.clf()
        plt.close()
        #reset parameters
        # sns.set_style(sns_default)
        ##save result to file
        self.save_result_tofile(rf_cv, X_test, y_test, cm_matrix, file_name=file_name,threshold=threshold_train)


        """
        #plt.show()
        #use the best parameters from cross validation redo RF
        rf = RandomForestClassifier(random_state=i,
                                    n_estimators=1000,max_features=0.5,class_weight="balanced",
                                    oob_score=True, n_jobs=8)
        rf.fit(X_train, y_train)
        
        y_pred = rf.predict(X_test)
        print("Train/test accuracy: {}/{}".format(
            accuracy_score(rf.predict(X_train), y_train),
            accuracy_score(rf.predict(X_test), y_test)))
        # roc for train data
        show_roc(self,rf, X_train, y_train, (file_name+"train"))
        #roc for test data
        show_roc(self,rf, X_test, y_test, (file_name + "test"))
        print("fprlength:{}".format(len(fpr)))

        # show the plot
        print("roc score:{}".format(roc_auc))
        rf_pro = rf.predict_proba(X_train)
        print(len(rf_pro)==len(y_train))
        neg_pro=[]
        pos_pro=[]
        y_train_pred=rf.predict(X_train)
        for i,line in enumerate(rf_pro):
            if y_train.iloc[i]==1:
                pos_pro.append(line[1])
            else:
                neg_pro.append(line[1])
            print("{}, label:{} prelabel:{}".format(line,y_train.iloc[i],y_train_pred[i]))
        print("mean probability for label 1 group: {}".format(sum(pos_pro)/len(pos_pro)))
        print("mean probability for predict label 0 as 1 group: {}".format(
            sum(neg_pro) / len(neg_pro)))
        print("min_for label1: {}".format(min(pos_pro)))
        print("max_for label2: {}".format(max(neg_pro)))

        plt.scatter(y=pos_pro,x=range(len(pos_pro)),c="red")
        plt.scatter(y=neg_pro, x=range(len(neg_pro)),c="green")
        plt.savefig("probability for two label{}".format(file_name))
        plt.show()
        threshold=(sum(pos_pro)/len(pos_pro))
        print(threshold)
        y_pre_threshold=[]
        for point in rf.predict_proba(X_test):
            if point[1]>threshold:
                y_pre_threshold.append(1)
            else:
                y_pre_threshold.append(0)
        else:
            cm = confusion_matrix(y_test, y_pre_threshold)
            cm_display = ConfusionMatrixDisplay(cm).plot()
            cm_display.figure_.savefig('cm_{}_{}.png'.format(threshold,file_name),
                                       dpi=300)
            plt.title("RF confusion matrix with best parameters")
            cm_display.figure_.savefig(
                'cm_threshold{}_{}.png'.format(threshold,file_name), dpi=300)
            #plt.show()

            plt.close()
        cm = confusion_matrix(y_test, y_pred)
        cm_display = ConfusionMatrixDisplay(cm).plot()
        cm_display.figure_.savefig('cm_{}.png'.format(file_name), dpi=300)
        plt.title("RF confusion matrix with best parameters")
        cm_display.figure_.savefig(
            'cm_cv_best_parameters{}.png'.format(file_name), dpi=300)
        #plt.show()

        plt.close()
        plt.figure()
        fi = pd.DataFrame(data=rf.feature_importances_, index=X_train.columns,
                          columns=['Importance']) \
            .sort_values(by=['Importance'], ascending=False)

        #And visualize
        sns.barplot(data=fi.head(20), x="Importance", y=(fi.head(20)).index).set_title(
            "feature importance for RF model")
        #fig = ax.get_figure()
        plt.savefig('../feature_importance{}.png'.format(file_name))
        plt.close()

        """
        # #save model
        filename = '../autodata/model/rf_test_model_cv{}'.format(file_name)
        joblib.dump(rf_cv, filename)
        return rf_cv

    def SVM(self,X_train, X_test, y_train, y_test,file_name="",i:int=0):
        # pars = [{'C': [ 0.1]},
        #         {'kernel':['poly','rbf','sigmoid']},
        #         {'class_weight':["balanced", None]}]
        # # The same probability calibration procedure is available for all estimators via the CalibratedClassifierCV (see Probability calibration).
        # # In the case of SVC and NuSVC, this procedure is builtin in libsvm which is used under the hood, so it does not rely on scikit-learn’s CalibratedClassifierCV.
        # svc = GridSearchCV(SVC(class_weight="balanced",random_state=0,probability=True),
        #                    pars, cv=5, scoring="accuracy",verbose=3,n_jobs=2)
        svc=SVC(cache_size=1024,kernel="linear",class_weight="balanced",random_state=0,probability=True)
        svc_result=svc.fit(X_train, y_train)
        print(svc_result)
        y_pred = svc.best_estimator_.predict(X_test)
        cm = confusion_matrix(y_test, y_pred)
        cm_display = ConfusionMatrixDisplay(cm).plot()
        #cm_display.figure_.savefig('cm_svc_best_parameters{}.png'.format(file_name), dpi=300)
        #plt.show()
        plt.close()
        plt.figure()
        print("svc.cv_results_:{}".format(svc.cv_results_))
        sns.lineplot(y=svc.cv_results_["mean_test_score"],
                     x=svc.cv_results_['param_C'].data,
                     hue=svc.cv_results_['param_kernel'])
        plt.xlabel("max features")
        plt.ylabel("Accuracy (mean 3-fold CV)")
        plt.title(
            "Accuracy with different estimators and features for SVC model")
        #plt.savefig("../Accuracy with different estimators and features for SVC model_{}".format(file_name))

    def save_result_tofile(self,rf,X_test, y_test,cm_matrix,file_name="",threshold=0):
        """

        :param rf:
        :param X_test:
        :param y_test:
        :param cm_matrix:
        :param file_name:
        :return:
        """
        y_pred = rf.best_estimator_.predict(X_test)
        y_probs = rf.predict_proba(X_test)
        # keep probabilities for the positive outcome only
        yhat = y_probs[:, 1]
        result=pd.DataFrame()
        print(yhat)
        result["yhat"] = yhat
        result["y_test"] = y_test
        result.to_csv('../autodata/model/{}_predictresult.csv'.format(file_name))
        ##distribution for the probability of prediction
        label_1_pre=[]
        label_0_pre=[]
        for i,index in enumerate(y_test.index):
            #print(y_test.loc[index])
            if y_test.loc[index]==1:
                label_1_pre.append(i)
            else:
                label_0_pre.append(i)

        self.distribution_plot(yhat[label_1_pre],yhat[label_0_pre],file_name,"prediction")


        accuracy_score_function = accuracy_score(y_test, y_pred)
        mcc_score = matthews_corrcoef(y_test, y_pred)
        roc_score = roc_auc_score(y_true=y_test,y_score=yhat)

        TN = cm_matrix[0][0]
        FP = cm_matrix[0][1]
        FN = cm_matrix[1][0]
        TP = cm_matrix[1][1]
        print(rf.best_score_)
        accuracy = ((TP + TN) / (TP + TN + FP + FN))
        print("TP:{}".format(TP))
        print("FP:{}".format(FP))
        print("FN:{}".format(FN))
        print("TN:{}".format(TN))
        print("accuracy:{}".format(((TP+TN)/(TP+TN+FP+FN))))
        sensitivity = TP / (TP + FN)
        specificity = TN / (TN + FP)
        precision = TP / (TP + FP)
        file1 = open(
            "../autodata/separate_seed_result/{}prediction_summary.txt".format(
                file_name), "a")
        file1.write("The following is the result for {}:".format(file_name))
        file1.write("best parameters:{}\n".format(rf.best_params_))
        file1.write("accuracy_score: {}\n".format(accuracy_score_function))
        file1.write("accuracy: {}\n".format(accuracy))
        file1.write("sensitivity : {}\n".format(sensitivity))
        file1.write("specificity : {}\n".format(specificity))
        file1.write("precision : {}\n".format(precision))
        file1.write("mcc_score: {}\n".format(mcc_score))
        file1.write("roc_score: {}\n".format(roc_score))
        file1.write("threshold: {}\n".format(threshold))
        #file1.write("best_score: {}\n".format(rf.best_score_))
        file1.write("###############################################")
        file1.close()

    def show_roc(self,rf,X,y,file_name):
        from sklearn import metrics

        # #save drfault parameters
        # print(plt.rcParams)
        # default_para=plt.rcParams
        # roc for training data
        y_probs = rf.predict_proba(X)
        # keep probabilities for the positive outcome only
        yhat = y_probs[:, 1]
        print(yhat)
        fpr, tpr, thresholds = metrics.roc_curve(y, yhat)
        roc_auc = metrics.auc(fpr, tpr)
        gmeans = []
        for point in range(len(fpr)):
            gmean = sqrt(tpr[point] * (1 - fpr[point]))
            gmeans.append(gmean)
        # # locate the index of the largest g-mean
        ix = np.argmax(gmeans)
        print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))

        display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr,
                                          roc_auc=roc_auc).plot()
        plt.title('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
        display.figure_.savefig("../autodata/separate_seed_result/RF_ROC_curve_{}_data.svg".format(file_name))
        #plt.show()
        plt.close()
        # #reset to rc parameter before
        # plt.rcParams.update(default_para)
        ###precision_recall###
        from sklearn.metrics import PrecisionRecallDisplay

        display = PrecisionRecallDisplay.from_estimator(
            rf, X, y, name="Random Forest"
        )
        precision_recall_figure = display.ax_.set_title("Precision-recall curve")
        display.figure_.savefig(
            "../autodata/separate_seed_result/RF_precision_recall_figure_{}_data.svg".format(
                file_name))
        #plt.show()


        plt.close()
        return thresholds[ix]

    def cm_threshold(self,threshold,x,y,rf,file_name):
        """
        The function is to plot confusion matrix with set threshold

        """
        #save drfault parameters
        print(plt.rcParams)
        default_para=plt.rcParams
        y_pre_threshold = []
        for point in rf.predict_proba(x):
            if point[1] >= threshold:
                y_pre_threshold.append(1)
            else:
                y_pre_threshold.append(0)
        else:
            #plt.rcParams["figure.figsize"] = (20, 20)
            #plt.rcParams.update({'font.size': 8})
            cm = confusion_matrix(y, y_pre_threshold)
            fig,ax=plt.subplots()
            #sns.set(font_scale=4)
            # sns.heatmap(cm, annot=True, annot_kws={"size": 16},ax=ax)  # font size
            # plt.xlabel("Predicted label")
            # plt.ylabel("True label")
            # plt.title(
            #     "RF confusion matrix threshold:{}".format(
            #         threshold))
            # plt.savefig("heatmap_cmplot{}.png".format(file_name))
            # plt.close()
            cm_display = ConfusionMatrixDisplay(cm).plot()
            # cm_display.figure_.savefig(
            #     'cm_{}_{}.png'.format(threshold, file_name),dpi=300)
            plt.title(
                "RF confusion matrix threshold:{}".format(
                    threshold))
            cm_display.figure_.savefig(
                '../autodata/separate_seed_result/cm_threshold{}_{}.svg'.format(
                    threshold, file_name), dpi=300)
            # plt.show()
            # plt.rcParams.update(default_para)
            plt.close()
        return cm

    def distribution_plot(self,label_1:list,label_0:list,file_name="",titile=""):
        """
        This is the function for probability distribution

        """
        # Draw the density plot
        #fig, ax = plt.subplots()
        print(plt.rcParams)
        default_para=plt.rcParams
        #se the same bin
        fig, ax =plt.subplots(figsize=(50,30))
        bin=[x*0.01 for x in list(range(1,101,1))]
        sns.histplot(label_1, kde=True,
                     color="orange",
                     label="label1",ax=ax,bins=bin)
        sns.histplot(label_0, kde=True,
                     color="blue",
                     label="label0",ax=ax,bins=bin)
        # Plot formatting
        plt.tight_layout()
        plt.legend(prop={'size': 8}, title='Label')
        plt.title('{}_{} probability distribution'.format(file_name,titile))
        plt.xlabel('probability')
        plt.savefig("../autodata/separate_seed_result/{} probability distribution.svg".format(file_name),bbox_inches = 'tight')
        #reset to the parameter before drawing the plot
        #plt.rcParams.update(default_para)
        plt.close()
        # plt.ylabel('Density')
        #plt.show()

    def hierarchical_clustering(self,sequence_data,label):

        from sklearn.cluster import AgglomerativeClustering
        from scipy.cluster.hierarchy import fcluster, cut_tree, linkage, \
            dendrogram
        import matplotlib.pyplot as plt


        # pca_df = self.run_PCA(sequence_data,
        #                        label, file_name="kmer_sequences")
        fig, axes = plt.subplots(2, 3, figsize=(15, 15))
        axis = [[0, 0], [0, 1], [0, 2],
                [1, 0], [1, 1], [1, 2]]
        methods = ['complete', 'average', 'single'] * 2
        metric = ['correlation', 'euclidean'] * 3
        # for metd, metr, j in zip(methods, metric, axis):
        #     # print(j[0], j[1], i)
        sequence_data.drop("methyl_type", inplace=True,axis=1)
        print(sequence_data.to_numpy())
        data = sequence_data.to_numpy()
        # hc = linkage(data, method=metd, metric=metr)
        # hc_clusters = cut_tree(hc, 4).ravel()
        cluster=AgglomerativeClustering(n_clusters=5,affinity="euclidean")
        cluster.fit(data)
        labels=cluster.fit_predict(data)
        print(labels)
        data["cluster_label"]=list(labels)
        data.to_csv("clustering_label.csv")
        # axes[0,0].scatter(pca_df.PC1, pca_df.PC2, c=hc_clusters,
        #                          s=5)
        # axes[0,0].set_title(f'Method = "euclidean"')
        # axes[0,0].set_xlabel('PC1')
        # axes[0,0].set_ylabel('PC2')
        # fig.suptitle('Hierarchical clusters on PCA 2 Dimension ')
        # plt.savefig("clustering")
        #plt.show()

    def predict(self,substrate_smile_df,enzyme_sequences:str="", inputmodel = None,num_bits: int = 2048):
        """

        :param substrate_smile_df:
        :param enzyme_sequences:
        :param inputmodel:
        :return:
        """
        if inputmodel:
            model = inputmodel
        else:
            if enzyme_sequences=="":
                #load model
                model = joblib.load('../data/model/rf_test_model')
        #create new columns
        substrate_smile_df["substrate_mol"] = pd.DataFrame(
            len(substrate_smile_df.index) * [0]).astype('object')
        for index in substrate_smile_df.index:
            substrate_smile = substrate_smile_df.loc[index,"substrate_smiles"]
            #smile to fingerprint
            try:
                substrate = Chem.MolFromSmiles(substrate_smile)
                #reset model index and mapnumber
                self_defined_mol_object = Molecule()
                substrate, no_use_virable = self_defined_mol_object.mol_with_atom_index(mol_object =substrate)

                substrate_smile_df.loc[index,"substrate_mol"] = substrate
            except:
                print("unable to proces the smile")
                print(substrate_smile)
                substrate_smile_df.loc[index, "substrate_mol"] = "NA"
        substrate_df = substrate_smile_df
        methyl_site_atom = {}

        #need to chang to prepare data based on model
        fingerprint_df=self.fingerprint_df_preapare(substrate_df,num_bits=num_bits,radius=3)
        Y = model.predict(fingerprint_df[list(fingerprint_df.columns)[:-2]])
        print(Y)

        for index,y in enumerate(Y):
            mol_id = fingerprint_df.loc[index,'molecular_id']
            if mol_id not in methyl_site_atom:
                methyl_site_atom[mol_id]=[]
            else:
                if y == 1:
                    atom_index=int(fingerprint_df.loc[index,'atom_index'])
                    atom = substrate.GetAtomWithIdx(atom_index)
                    methyl_site_atom[mol_id].append((atom.GetSymbol()+str(atom_index)))
                else:
                    continue
        else:
            return methyl_site_atom

    def duplicate_1_class(self,input,times):

        indexNames = input[input['label'] == 1].index
        methyl_site_data = input.loc[indexNames]
        for i in range(times):
            input = input.append(methyl_site_data)
        print(input)
        return copy.deepcopy(input)

    def multidata_predict(self,dataset=None,inputmodel = None):

        if (dataset == None) or (inputmodel == None):
        #then do with test data
            loded_data=pd.read_excel("../data/prediction_test.xlsx", header=0, index_col=0)
            print(loded_data)
            methyl_site_atom= self.predict(loded_data,num_bits=1024)
            print(methyl_site_atom)

    def use_less_similar_data_for_test(self,df,test:float=0.2,group_column:str='molecular_id',i:int=0,num_bit=2048):
        """

        :param df: input_dataframe with substrate fingerprint and sequences encoding
        :param test: the propotion of data use for test data
        :param group_column: group based on the column
        :param i: int, random_state
        :param num_bit: substrate fingerprint length
        :return:
        """

        splitter = GroupShuffleSplit(test_size=test, n_splits=1,
                                     random_state=i)
        split = splitter.split(df, groups=df[group_column])
        train_inds, test_inds = next(split)
        train = df.iloc[train_inds]
        test = df.iloc[test_inds]

        simliar_dictionary=self.create_similarity_dictionary()
        self.check_test_train_similarity(test,
                                         train,
                                         num_bit,similarity_dictionary=simliar_dictionary)

        train=self.train
        print('Q6AWU6' in list(train["Entry"]))
        print('Q6AWU6' in list(test["Entry"]))
        X_train = (copy.deepcopy(train)).drop(
            columns=["Entry", "molecular_id", "label"])
        X_test = (copy.deepcopy(test)).drop(
            columns=["Entry", "molecular_id", "label"])
        Y_train = train["label"]
        Y_test = test["label"]


        # print(X_train, X_test, Y_train, Y_test)

        return X_train, X_test, Y_train, Y_test

    def check_test_train_similarity(self, test, train, numbit,similarity_dictionary:dict={},maxmum_similarity:float=1.0,minimum_similarity:float=0.1):
        """

        :param test: entry in test
        :param train: entry in train
        similarity_dictionary: {entry；{entry1:identity1,entry2:identity2...}}
        :return: train with sequences within the range of similarity
        """

        #print(len(test.index))
        # this dictionary is to save fingerprint and corresponding index in train data
        #print(train.index)
        train_sub_fg_dictionary={}
        for index_train in train.index:
            column_name = list(train.columns)[:numbit * 2]
            #column_name=[x for x in range(numbit * 2)]
            list_sub_fg = (train.loc[index_train,column_name]).values.tolist()
            #print(list_sub_fg)
            list_sub_fg = [str(x) for x in list_sub_fg]
            sub_finderprint = "".join(list_sub_fg)
            #print(sub_finderprint)
            if sub_finderprint not in train_sub_fg_dictionary.keys():
                train_sub_fg_dictionary[sub_finderprint]=[index_train]
            else:
                train_sub_fg_dictionary[sub_finderprint].append(index_train)

        #list odf all fingerprint present in train data
        train_sub_fg_list=list(train_sub_fg_dictionary.keys())

        keep_index=[]
        for index_test in test.index:
            #join the fingerprint to a sting for searching
            column_name = list(test.columns)[:numbit * 2]
            list_sub_fg = (test.loc[index_test,column_name]).values.tolist()
            list_sub_fg = [str(x) for x in list_sub_fg]
            sub_fingerprint = "".join(list_sub_fg)
            if sub_fingerprint in train_sub_fg_list:
                test_entry = test.loc[index_test,"Entry"]
                #test_seq=seq_dictionary[test_entry]
                #check all the index for this substrate in train data
                for index_train in train_sub_fg_dictionary[sub_fingerprint]:
                    #print(test_entry,"test",test_entry)
                    train_entry = train.loc[index_train,"Entry"]
                    #train_seq=seq_dictionary[train_entry]
                    #if test seq is similar to test seq, remove this row in train
                    if train_entry in similarity_dictionary.keys():
                        if test_entry in similarity_dictionary[train_entry].keys():
                            #check the similarity
                            similarity = similarity_dictionary[train_entry][test_entry]

                            if (similarity < maxmum_similarity) and (similarity > minimum_similarity):
                                # print(similarity_dictionary[train_entry][
                                #           test_entry])
                                # print(maxmum_similarity)
                                # print(minimum_similarity)
                                keep_index.append(index_test)
                    else:
                        if minimum_similarity==0.0:
                            keep_index.append(index_test)
                        #print("There are no similar sequences with {}".format(test_entry))
        #print(keep_index)
        keep_index = list(set(keep_index))
        test= test.loc[keep_index]
        #test = pd.DataFrame(test,index=keep_index)
        # print(len(test.index))
        # print(test)
        # self.train = train
        return test

    def create_similarity_dictionary(self,file="../autodata/all_seq_smilarity_result.tab"):
        """
        This is the function to get the similarity dictionary

        :param file: string, file name and path with mmseqs result
        :return:
        similarity_dictionary: {entry1；{entry2:identity2,...},entry2:{...}}
        """

        pd_result_mmseqs=pd.read_csv(file,header=None,index_col=None,sep="\t")
        pd_result_mmseqs=pd_result_mmseqs.loc[:,0:2]
        pd_result_mmseqs.columns=["Entry1","Entry2","identity"]
        pd_result_mmseqs=pd_result_mmseqs.sort_values('identity', ascending=False)

        # remove sequences compared to themselves and seed sequences
        for index in pd_result_mmseqs.index:
            if pd_result_mmseqs.loc[index,"Entry1"] == pd_result_mmseqs.loc[index,"Entry2"]:
                pd_result_mmseqs.drop(index=index,inplace=True)
            elif "-" in (pd_result_mmseqs.loc[index,"Entry1"] or pd_result_mmseqs.loc[index,"Entry2"]):
                pd_result_mmseqs.drop(index=index,inplace=True)

        # print(pd_result_mmseqs)
        # print(len(pd_result_mmseqs.index))
        count = (pd_result_mmseqs['identity'] > 0.8).sum()
        print(count)
        similarity_dictionary = {}
        for i in pd_result_mmseqs.index:
            entry1 = pd_result_mmseqs.loc[i,"Entry1"]
            entry2 = pd_result_mmseqs.loc[i,"Entry2"]
            identity = pd_result_mmseqs.loc[i,"identity"]
            # save in the dictionary idf identity larger than 0.5
            if entry1 not in similarity_dictionary.keys():
                similarity_dictionary[entry1]={entry2:identity}
            else:
                similarity_dictionary[entry1][entry2] = identity


        print(similarity_dictionary)
        #save similarity dictionary for later use
        with open("../autodata/similarity_dictionary", "wb") as dill_file:
            dill.dump(similarity_dictionary, dill_file)

        return similarity_dictionary

    def evaluate_model_performance_based_on_molecule(self,input_file):
        #load model
        model = joblib.load('../data/model/rf_test_model')
        #load test data

        input_dataframe = pd.read_csv(input_file,header=0, index_col=0)
        print(input_dataframe)
        X_train, X_test, y_train, y_test = self.prepare_train_teat_data(
            input_dataframe)
        y_predict=model.predict(X_test)

    def comapre_result_for_different_similarity_test(self,input_file):

        #load similarity dictionary
        try:
            with open("../autodata/similarity_dictionary", 'rb') as file1:
                similarity_dictionary = dill.load(file1)
        except:
            print("create similarity dictionary")
            similarity_dictionary = self.create_similarity_dictionary()
        #print(similarity_dictionary)
        #split test data bsed on different sequences similarity
        input_df=pd.read_csv(input_file,header=0,index_col=0)
        splitter = GroupShuffleSplit(test_size=0.4, n_splits=1, random_state=0)
        split = splitter.split(input_df, groups=input_df["molecular_id"])
        train_inds, test_inds = next(split)
        for similarity in [(0.0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8),(0.8,1.0)]:
            train = input_df.iloc[train_inds]
            test = input_df.iloc[test_inds]
            print("test len {}".format(len(test)))
            print(similarity)
            minmum,maxmum = similarity
            test=self.check_test_train_similarity(test, train, 166,similarity_dictionary=similarity_dictionary,minimum_similarity=minmum,maxmum_similarity=maxmum)
            print(len(test))
            X_train = (copy.deepcopy(train)).drop(columns=["Entry","label", "atom_index"])
            Y_train = train["label"]
            X_test = (copy.deepcopy(test)).drop(columns=["Entry", "label", "atom_index"])
            X_test.to_csv("167fg_bond3_rf{}_{}X_test.csv".format('MACCS',similarity))
            Y_test = test["label"]
            #print(Y_test.index)
            Y_test.to_csv("167fg_bond3_rf{}_{}Y_test.csv".format('MACCS',similarity))
            X_train=X_train.drop(columns=["molecular_id","methyl_type"])
            X_test=X_test.drop(columns=["molecular_id","methyl_type"])
            #print(X_test)
            self.RF_model(X_train, X_test, y_train=Y_train, y_test=Y_test,
                                file_name="167fg_bond3_rf{}_{}%_{}%".format('MACCS',str(int(minmum*100)),str(int(maxmum*100))),i=0)

def create_MACCSkey_fingerprint():
    model=Model_class()
    model.check_file_exist("../autodata/seq_smiles_all.csv","../autodata/diction_atom_all")
    data_with_site = pd.read_csv("../autodata/seq_smiles_all.csv",
                                 header=0, index_col=0)
    with open('../autodata/diction_atom_all', 'rb') as file1:
        diction_atom = dill.load(file1)
    #model.save_MACCSkeys_fingerprints_to_dataframe(data_with_site, diction_atom,3, drop_atoms=True,file_name="all_data")
    model.save_fingerprints_to_dataframe(data_with_site, diction_atom,num_bits=2048,radius=3, drop_atoms=True,file_name="all_data")

def main():
    #create_MACCSkey_fingerprint()

    model = Model_class()
    model.comapre_result_for_different_similarity_test("../autodata/input_data/active_site/PF08241_bit_score11_coverage0.5_ACS_bit167_3_remove_redundant_MACCS.csv")
if __name__ == "__main__":
    main()
