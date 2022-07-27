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
from script.molecular_class import molecular, reaction
from script import parse_data
#for data structure
import pandas as pd
import numpy as np

from sys import argv
from urllib.request import urlopen
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint, SDMolSupplier
from rdkit.Chem.Draw import IPythonConsole
import pikachu
from rdkit.Chem import Draw, AllChem, rdmolops
import glob
import joblib
#for display and draw picture and graph
from IPython.display import SVG, display, Image
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

from PIL import Image

#import for model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import confusion_matrix, mean_squared_error, mean_absolute_percentage_error, r2_score, mean_absolute_error, accuracy_score, ConfusionMatrixDisplay, multilabel_confusion_matrix
from sklearn.model_selection import GridSearchCV
def merge_uniprot_id_smile(rheauniprot_dataframe,seq_df):
    """
    combine sequence and substrate smile in onefile and save to csv file
    :param
        rheauniprot_dataframe:
        seq_df:
    :return:
    """
    df1 = parse_data.get_substrate_chebi("data/Rhea-ec_2_1_1.tsv")
    df1 = parse_data.get_smile(df1)
    df1.index = df1.index.map(lambda x: x.split(":")[1])
    df1["RHEA_ID"] = df1.index
    df1.set_index("RHEA_ID", inplace=True)
    uniprot_entry = parse_data.remove_duplicated_id(
        r"E:\Download\regioselectivity_prediction\data\hmm_out")
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


def drop_useless_column(dataframe):
    dataframe = dataframe.loc[:,["Sequence","sub_mols","pro_mols","sub_smiles","pro_smiles"]]
    #print(dataframe)


def keep_longest_smile(dataframe_before):

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


def return_reactions(dataframe_rr):
    """

    :param dataframe:
    :return:
    """
    atom_object_dictionary = {}
    dataframe_rr["reactant_site"] = pd.DataFrame(
        len(dataframe_rr.index) * [0]).astype('object')
    dataframe_rr["mainsub_mol"] = pd.DataFrame(
        len(dataframe_rr.index) * [0]).astype('object')
    for index in dataframe_rr.index:
        rxn = dataframe_rr.loc[index,"rxn"]
        sub = dataframe_rr.loc[index, "main_sub"]
        pro = dataframe_rr.loc[index, "main_pro"]
        reaction1 = reaction(substrates=sub, products=pro)
        #print(dataframe.loc[index,"RHEA_ID"])
        rxn_file_name = "data/rxn_picture/{}".format(dataframe_rr.loc[index,"Entry"])
        #r1 = reaction1.get_reaction_sites(rxn_object=rxn,file_name=rxn_file_name)
        r2, index_list, mainsub_mol = reaction1.get_reactant_atom()
        print("index: {}".format(index))

        #save atom index(methylation site) from substrate in dataframe
        site = ",".join(index_list)
        print(site)
        dataframe_rr.loc[index,"mainsub_mol"] = mainsub_mol
        if site != "":
            dataframe_rr.loc[index, "reactant_site"] = site
        else:
            #because we use the largest molecular, but for some substrates, methyl donor is larger
            #we first leave those
            dataframe_rr.loc[index,"reactant_site"] = "NA"
        #save list of methyl site atom objects and index to a dictionary
        atom_object_dictionary[index] = site
    else:
        with open("data/seq_smiles", "wb") as dill_file:
            dill.dump(dataframe_rr, dill_file)
        with open("data/diction_atom", "wb") as dill_file:
            dill.dump(atom_object_dictionary, dill_file)

        return dataframe_rr,atom_object_dictionary

def save_fingerprints_to_dataframe(sauce_data,atom_object_dictionary,num_bits: int = 2048,radius: int = 3):
    """
    this function is to build inputdata with fingerprints and labels
    :param sauce_data:
    :param atom_object_dictionary:
    :param num_bits:
    :param radius:
    :return:
    """
    self_defined_mol_object = molecular()
    input_dataframe = pd.DataFrame()
    current_index = 0
    print(sauce_data)
    for index in sauce_data.index:
        print(index)
        sub_mol = sauce_data.loc[index,"mainsub_mol"]
        #Draw.ShowMol(sub_mol, size=(600, 600))
        sub_rest_mol, no_use_variable = self_defined_mol_object.mol_with_atom_index(mol_object=sub_mol)
        fingerprint_mol = self_defined_mol_object.create_fingerprint_mol(
            sub_rest_mol, num_bits=num_bits, radius=radius)
        for atom in sub_mol.GetAtoms():
            #set label
            sy_index = (atom.GetSymbol() + str(atom.GetIdx())+":"+str(atom.GetAtomMapNum()))
            if sy_index in atom_object_dictionary[index]:
                label = 1
            else:
                label = 0
            #atom_index_sub = atom.GetAtomMapNum()
            newrow = {}
            #resrt atom index and then build fingerprint
            fingerprint_atom = self_defined_mol_object.create_fingerprint_atom(
                sub_rest_mol, atom_object=atom, num_bits=num_bits, radius=radius)
            #add fingerprint atom index ebedding sequences and label to dataframe
            for i,item in enumerate(fingerprint_mol):
                newrow[i] = item
            for j,item in enumerate(fingerprint_atom):
                newrow[j+i] = item
            #newrow['atom_index'] = atom_index_sub
            add_dataframe = pd.DataFrame(newrow,index=[current_index])
            input_dataframe=pd.concat([input_dataframe,add_dataframe],axis=
                                      0)
            input_dataframe.loc[current_index, "molecular_id"] = "m"+str(index)
            input_dataframe.loc[current_index,"label"] = label
            current_index += 1
    print(input_dataframe)
    input_dataframe.to_csv("data/input_dataframe_withoutstructure_{}.csv".format(str(num_bits)))
    with open("data/input_dataframe_withoutstructure_{}".format(str(num_bits)), "wb") as dill_file:
        dill.dump(input_dataframe, dill_file)

    return input_dataframe
#
# def sequences_feature_to_dataframe(substratedata,whole_data):
#     for index in substratedata.index:
#

def fingerprint_df_preapare(substrate_mol_df,num_bits: int = 2048,radius: int = 3):
    """

    :param substrate_mol:
    :param num_bits:
    :param radius:
    :return:
    """
    current_index = 0
    mol_id:int = 0
    fingerprint_dataframe = pd.DataFrame()
    self_defined_mol_object = molecular()
    print(substrate_mol_df)
    for index in substrate_mol_df.index:
        substrate_mol = substrate_mol_df.loc[index,'substrate_mol']
        try:

            #Draw.ShowMol(substrate_mol, size=(600, 600))
            fingerprint_mol = self_defined_mol_object.create_fingerprint_mol(
                substrate_mol, num_bits=num_bits, radius=radius)
            for atom in substrate_mol.GetAtoms():
                fingerprint_atom = self_defined_mol_object.create_fingerprint_atom(
                    substrate_mol, atom_object=atom, num_bits=num_bits, radius=radius)
                atom_index_sub = atom.GetIdx()
                newrow = {}
                # add fingerprint atom index ebedding sequences and label to dataframe
                for i, item in enumerate(fingerprint_mol):
                    newrow[i] = item
                for j, item in enumerate(fingerprint_atom):
                    newrow[j + i] = item

                add_dataframe = pd.DataFrame(newrow, index=[current_index])
                fingerprint_dataframe = pd.concat([fingerprint_dataframe, add_dataframe], axis=
                0)
                fingerprint_dataframe.loc[current_index, "molecular_id"] = "m"+str(mol_id)
                fingerprint_dataframe.loc[current_index, 'atom_index'] = atom_index_sub
                current_index += 1
        except:
            print("skip in mol{}".format(mol_id))
        mol_id +=1
    print(fingerprint_dataframe)
    return fingerprint_dataframe


def prepare_train_teat_data(df,test:float=0.2,group_column:str='molecular_id'):
    """
    simply sepreate train and test data
    :param data:
    :return:
    """

    from sklearn.model_selection import GroupShuffleSplit
    splitter = GroupShuffleSplit(test_size=test, n_splits=1, random_state=1)
    split = splitter.split(df, groups=df[group_column])
    train_inds, test_inds = next(split)
    train = df.iloc[train_inds]
    test = df.iloc[test_inds]

    X_train = train[list(train.columns)[:-2]]
    Y_train = train["label"]
    X_test = test[list(test.columns)[:-2]]
    Y_test = test["label"]

    # X = data[list(data.columns)[:-1]]
    # y = data["label"]
    #
    # train_test_split(X, y, random_state=1, stratify=y)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1,stratify=y)
    return X_train, X_test, Y_train, Y_test

def RF_model(X_train, X_test, y_train, y_test):
    """

    :param X_train:
    :param X_test:
    :param y_train:
    :param y_test:
    :return: randomforest model
    """
    # hyperparameters = {'n_estimators': [300, 30],
    #                    'max_features': [0.3,0.5,1.0]
    #                    }
    #
    # rf_cv = GridSearchCV(RandomForestClassifier(random_state=0),
    #                      hyperparameters,
    #                      cv=10,
    #                      verbose=True,
    #                      n_jobs=-1)
    #
    # rf_cv.fit(X_train, y_train)
    # print(rf_cv.best_params_)
    #
    # y_pred = rf_cv.best_estimator_.predict(X_test)
    #

    #use the best parameters from cross validation redo RF
    rf = RandomForestClassifier(random_state=1,
                                n_estimators=30,max_features=0.3,
                                oob_score=True)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)
    print("Train/test accuracy: {}/{}".format(
        accuracy_score(rf.predict(X_train), y_train),
        accuracy_score(rf.predict(X_test), y_test)))
    print()
    cm = confusion_matrix(y_test, y_pred)
    cm_display = ConfusionMatrixDisplay(cm).plot()
    cm_display.figure_.savefig('cm_128bit.png', dpi=300)
    plt.title("RF confusion matrix with best parameters")
    plt.show()

    fi = pd.DataFrame(data=rf.feature_importances_, index=X_train.columns,
                      columns=['Importance']) \
        .sort_values(by=['Importance'], ascending=False)

    # And visualize
    plt.figure(figsize=(200, 600))
    ax = sns.barplot(data=fi, x="Importance", y=fi.index).set_title(
        "feature importance for RF model")
    fig = ax.get_figure()
    plt.show()
    fig.savefig('feature_importance.png')
    #save model
    filename = 'data/model/rf_test_model'
    joblib.dump(rf, filename)
    return rf


def predict(substrate_smile_df,enzyme_sequences:str="", inputmodel = None):

    if inputmodel:
        model = inputmodel
    elif enzyme_sequences=="":
        #load model
        model = joblib.load('data/model/rf_test_model')
    #create new columns
    substrate_smile_df["substrate_mol"] = pd.DataFrame(
        len(substrate_smile_df.index) * [0]).astype('object')
    for index in substrate_smile_df.index:
        substrate_smile = substrate_smile_df.loc[index,"substrate_smiles"]
        #smile to fingerprint
        try:
            substrate = Chem.MolFromSmiles(substrate_smile)
            #reset model index and mapnumber
            self_defined_mol_object = molecular()
            substrate, no_use_virable = self_defined_mol_object.mol_with_atom_index(mol_object =substrate)

            substrate_smile_df.loc[index,"substrate_mol"] = substrate
        except:
            print("unable to proces the smile")
            print(substrate_smile)
            substrate_smile_df.loc[index, "substrate_mol"] = "NA"
    substrate_df = substrate_smile_df
    methyl_site_atom = {}
    #need to chang to prepare data based on model
    fingerprint_df=fingerprint_df_preapare(substrate_df,128,3)
    Y = model.predict(fingerprint_df[list(fingerprint_df.columns)[:-2]])
    print(Y)

    for index,y in enumerate(Y):
        mol_id = fingerprint_df.loc[index,'molecular_id']
        if mol_id not in methyl_site_atom:
            methyl_site_atom[mol_id]=[]
        else:
            if y == 1:
                atom_index=fingerprint_df.loc[index,'atom_index']
                atom = substrate.GetAtomWithIdx(atom_index)
                methyl_site_atom[mol_id].append((atom.GetSymbol()+str(atom_index)))
            else:
                continue
    else:
        return methyl_site_atom

def multidata_predict(dataset=None):
    if dataset == None:
    #then do with test data
        loded_data=pd.read_excel("data/prediction_test.xlsx", header=0, index_col=0)
        print(loded_data)
        methyl_site_atom= predict(loded_data)
        print(methyl_site_atom)
def main():
    # readfile whihc contains the rhea id and related uniprotid
    #run with command line
    # rh_file = argv[1]
    #run with pycharm
    '''
    rh_file = "data/rhea2uniprot_sprot.tsv"
    rheauniprot_dataframe = parse_data.readrhlist(rh_file)

    #read id and sequences in dataframe
    seq_file = "data/id_tosequence.xlsx" # run with pycharm
    #seq_file = argv[2]
    id_seq_dataframe = parse_data.read_sequence(seq_file)
    seq_smiles = merge_uniprot_id_smile(rheauniprot_dataframe,id_seq_dataframe)
    data_frame = keep_longest_smile(seq_smiles)
    data_with_site,diction_atom = return_reactions(data_frame)
    print(data_with_site["reactant_site"])
    '''
    with open('data/seq_smiles','rb') as file1:
        data_with_site = dill.load(file1)
    with open('data/diction_atom','rb') as file1:
        diction_atom = dill.load(file1)
    indexNames = data_with_site[data_with_site['reactant_site'] == 'NA'].index
    # Delete these row indexes from dataFrame
    data_with_site.drop(indexNames, inplace=True)
    print(len(data_with_site.index))
    save_fingerprints_to_dataframe(data_with_site,diction_atom,56,3)

    #read manual_data
    #parse_data.read_mannual_data()







    # X = pd.read_csv("data/input_dataframe_withoutstructure_128.csv",header=0,index_col=0)
    # X_train, X_test, y_train, y_test = prepare_train_teat_data(X)
    # #train RF model
    # model = RF_model(X_train, X_test, y_train, y_test)
    multidata_predict()
    #print(id_seq_dataframe)
    #link the sequences and reaction participant put in one csv file
    #fulldata = id_seq_dataframe.join(rheauniprot_dataframe)
    # print(fulldata)
    # groups_uniprotid = fulldata.groupby(["Entry"])
    # for group in groups_uniprotid:
    #     print(group)
    #save to file
    # file_name = input("please input file name and directory")
    # fulldata.to_csv(file_name)
    #convert rheaid to uniprot entry, as there are no direct link between Rhea id to sequences
    #read rhea id file
    #rhea = open("data/Rhea-ec_2.1.1.list").readlines()
    #uniprot_list = rheaid_touniprot(rhea, rheauniprot_dataframe)
    #print(uniprot_list)
    # file = open("uniprot_list.txt","w")
    # for index,entry in enumerate(uniprot_list):
    #     if isinstance(entry, str):
    #         file.write("{},{}\n".format(entry,rhea[index]))
    #     else:
    #         entry.to_csv("uniprot_list.txt", mode='a',header=False)
    # print(uniprot_list)
    '''
    read_object_from_file()
    '''
    #if the nmbering from rdkit is following the carbon numbering rules?
    # molclass = molecular()
    # smile1='C1CC2=NC1=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=C2)N5'
    # smile2="C1=CC(=C(C=C1C=CC(=O)O)O)O"
    # mol,index = molclass.mol_with_atom_index(smile=smile2)
    # smile_withindex = Chem.MolToSmiles(mol)
    # for atom in mol.GetAtoms():
    #     atom.SetIsotope(0)
    # Draw.ShowMol(mol,size=(600,600))
    # smiles_with_atom_mappings = Chem.MolToSmiles(mol)



    #remove_duplicated_id("data/hmm_out/top_ten_hits_exclude_nomethylrelated/hmmscantbout_top_6_hit.pfam")

    # uniprot_entry = remove_duplicated_id(r"E:\Download\regioselectivity_prediction\data\hmm_out")
    # rheaid_to_uniprot(uniprot_entry, rheauniprot_dataframe)
    # draw_smile = pikachu.general.draw_smiles(smile_withindex)

    # rxn = AllChem.ReactionFromSmarts(r"[C:1]-[C:2]>>[C:1].[C:2]")
    # img = Draw.ReactionToImage(rxn,returnPNG=True)
    # #?
    # with open("test.png",'wb') as handle:
    #     handle.write(img)


if __name__ == "__main__":
    main()
