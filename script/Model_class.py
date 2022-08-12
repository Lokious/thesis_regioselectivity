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


class Model_class():

    #def __init__(self):
        #self.check_file_exist()

    def check_file_exist(self):

        try:
            file1 = 'data/seq_smiles_all'
            file2 = 'data/diction_atom_all'

            assert pathlib.Path(file1).exists()==True
            assert pathlib.Path(file2).exists()==True

        except:
            print("missing file{},{}".format(file1,file2))
            print("preparing the file......")
            rh_file = "data/rawdata/rhea2uniprot_sprot.tsv"
            rheauniprot_dataframe = parse_data.readrhlist(rh_file)
            #read id and sequences in dataframe
            seq_file = "data/rawdata/id_tosequence.xlsx" # run with pycharm
            id_seq_dataframe = parse_data.read_sequence(seq_file)
            seq_smiles = self.merge_uniprot_id_smile(rheauniprot_dataframe,id_seq_dataframe)
            #data_frame = self.keep_longest_smile(seq_smiles)
            data_frame=self.keep_methyled_substrate(seq_smiles)
            data_frame.to_csv("data/seq_dataframe_new.csv")
            self.return_reactions(data_frame)



    def merge_uniprot_id_smile(self,rheauniprot_dataframe,seq_df):
        """
        combine sequence and substrate smile in onefile and save to csv file
        :param
            rheauniprot_dataframe:
            seq_df:
        :return:
        """
        df1 = parse_data.get_substrate_chebi("data/rawdata/Rhea-ec_2_1_1.tsv")
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


    def drop_useless_column(self,dataframe):

        dataframe = dataframe.loc[:,["Sequence","sub_mols","pro_mols","sub_smiles","pro_smiles"]]
        #print(dataframe)


    def keep_longest_smile(self,dataframe_before):

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
        exclude H+,
        :param dataframe_before:
        :return:
        """
        #UNFINISHED
        dataframe_before["main_sub"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        dataframe_before["main_pro"] = pd.DataFrame(
            len(dataframe_before.index) * [0])
        for index in dataframe_before.index:
            subs = (dataframe_before.loc[index,"sub_smiles"]).split(".")
            pros = (dataframe_before.loc[index,"pro_smiles"]).split(".")
            for sub in subs:
                mol = Chem.MolFromSmiles(sub)
                if len(mol.GetAtoms())<2:
                    subs.remove(sub)
            for pro in pros:
                mol = Chem.MolFromSmiles(pro)
                if len(mol.GetAtoms())<2:
                    pros.remove(pro)
            if len(subs) == len(pros):
                reaction_object = reaction()
                smiles = reaction_object.main_substrate(subs, pros)
                if smiles:
                    main_sub, main_pro = smiles[0],smiles[1]
                    dataframe_before.loc[index, "main_sub"] = main_sub
                    dataframe_before.loc[index, "main_pro"] = main_pro

        return copy.deepcopy(dataframe_before)

    def return_reactions(self,dataframe_rr):
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
            if (sub != pro) and sub !=0 and pro !=0:
                reaction1 = reaction(substrates=sub, products=pro)
                #print(dataframe.loc[index,"RHEA_ID"])
                #rxn_file_name = "data/rxn_picture/{}".format(dataframe_rr.loc[index,"Entry"])
                #r1 = reaction1.get_reaction_sites(rxn_object=rxn,file_name=rxn_file_name)
                print("index: {}".format(index))
                r2, index_list, mainsub_mol = reaction1.get_reactant_atom()

                Draw.MolToFile(mainsub_mol,"data/substrate_mol_picture/{}.png".format(index),size=(600,600))
                #save atom index(methylation site) from substrate in dataframe
                site = ",".join(index_list)
                print(site)
                dataframe_rr.loc[index,"mainsub_mol"] = mainsub_mol

            else:
                site = ""
            if site != "":
                dataframe_rr.loc[index, "reactant_site"] = site
            else:
                #because we use the largest molecular, but for some substrates, methyl donor is larger
                #we first leave those
                dataframe_rr.loc[index,"reactant_site"] = "NA"
            #save list of methyl site atom objects and index to a dictionary
            atom_object_dictionary[index] = site
        else:
            with open("data/seq_smiles_all", "wb") as dill_file:
                                dill.dump(dataframe_rr, dill_file)
            with open("data/diction_atom_all", "wb") as dill_file:
                dill.dump(atom_object_dictionary, dill_file)
            dataframe_rr.to_csv("data/diction_atom_all.csv")


    def save_fingerprints_to_dataframe(self,sauce_data,atom_object_dictionary,num_bits: int = 2048,radius: int = 3,drop_atoms=False,file_name=""):
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

        for index in sauce_data.index:

            print(index)
            sub_mol = sauce_data.loc[index,"mainsub_mol"]
            # print(sub_mol)
            #Draw.ShowMol(sub_mol, size=(600, 600))
            sub_rest_mol, no_use_variable = self_defined_mol_object.mol_with_atom_index(mol_object=copy.deepcopy(sub_mol))
            fingerprint_mol = self_defined_mol_object.create_fingerprint_mol(
                sub_rest_mol, num_bits=num_bits, radius=radius)
            for atom in sub_mol.GetAtoms():


                #set label

                sy_index =(atom.GetSymbol() + str(atom.GetAtomMapNum()) + ":" + str(atom.GetIsotope()))
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
                            newrow[j + i] = item
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
                        input_dataframe.loc[current_index, "methyl_type"] = \
                        sauce_data.loc[index, "methyl_type"]
                        current_index += 1
                else:
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
                    input_dataframe.loc[current_index,"Entry"] = sauce_data.loc[index,"Entry"]
                    input_dataframe.loc[current_index, "methyl_type"] = \
                        sauce_data.loc[index, "methyl_type"]
                    current_index += 1
        print(input_dataframe)
        if drop_atoms:
            input_dataframe.to_csv("data/input_dataframe_withoutstructure_dropatoms{}.csv".format(file_name))
            with open("data/input_dataframe_withoutstructure_dropatoms{}".format(file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        else:
            input_dataframe.to_csv("data/input_dataframe_withoutstructure_{}.csv".format(file_name))
            with open("data/input_dataframe_withoutstructure_{}".format(file_name), "wb") as dill_file:
                dill.dump(input_dataframe, dill_file)
        return input_dataframe


    def fingerprint_df_preapare(self,substrate_mol_df,num_bits: int = 2048,radius: int = 3):
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


    def prepare_train_teat_data(self,df,test:float=0.2,group_column:str='molecular_id',i:int=0):
        """
        simply sepreate train and test data
        :param data:
        :return:
        """


        splitter = GroupShuffleSplit(test_size=test, n_splits=1, random_state=i)
        split = splitter.split(df, groups=df[group_column])
        train_inds, test_inds = next(split)
        train = df.iloc[train_inds]
        test = df.iloc[test_inds]

        #X_train = train[list(train.columns)[:-2]]
        X_train = (copy.deepcopy(train)).drop(columns=["Entry","molecular_id","label"])
        Y_train = train["label"]
        #X_test = test[list(test.columns)[:-2]]
        X_test = (copy.deepcopy(test)).drop(columns=["Entry", "molecular_id", "label"])
        Y_test = test["label"]
        #print(X_train, X_test, Y_train, Y_test)
        return X_train, X_test, Y_train, Y_test



    def _get_colors(self,num_colors):
        import colorsys
        colors = []
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
        import matplotlib.patches as mpatches
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
            title='Sum of first 2 loading value of PCA')

        plt.savefig("Sum_of_first_2_loading_value_of_PCA_{}".format(file_name))

        plt.close()
        pca_df = pd.DataFrame(pca_fit.fit_transform(data_with_site),index=data_with_site.index, columns=PC)

        fig, ax = plt.subplots()
        plt.scatter(pca_df.PC1, pca_df.PC2, c=colour_label, s=5)

        #show explained variance
        handle_list = []
        for key in map_dictionary.keys():
            handle_list.append(mpatches.Patch(color=map_dictionary[key], label=key))
        ax.legend(
            handles=handle_list)

        #plt.scatter(pca_df.PC1, pca_df.PC2,  s=5,c=colour_label)
        plt.xlabel("pc1")
        plt.ylabel("pc2")
        plt.title("First two component of PCA coloured by methylation type")
        plt.savefig(
            "pca_for encoding sequences and fingerprint_{}".format(file_name))
        plt.clf()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 pca_fit.explained_variance_ratio_, '-ro')
        plt.ylabel('Proportion of Variance Explained_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            'Proportion of Variance Explained_{}'.format(file_name))

        plt.clf()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 np.cumsum(pca_fit.explained_variance_ratio_), '-o')
        plt.ylabel('Cumulative Proportion of Variance Explained_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            'Cumulative Proportion of Variance Explained_{}'.format(file_name))

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
        plt.xlabel("pc1")
        plt.ylabel("pc2")
        plt.title("First two component of PCA coloured by label type")
        plt.savefig(
            "pca_for encoding sequences and fingerprint for label_{}".format(file_name))
        #plt.show()
        plt.close()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 pca_fit.explained_variance_ratio_, '-ro')
        plt.ylabel('Proportion of Variance Explained for label_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            'Proportion of Variance Explained_{}'.format(file_name))
        #plt.show()
        plt.clf()
        plt.plot(list(range(1, len(pca_df.columns) + 1)),
                 np.cumsum(pca_fit.explained_variance_ratio_), '-o')
        plt.ylabel(
            'Cumulative Proportion of Variance Explained foe label_{}'.format(file_name))
        plt.xlabel("components")
        plt.savefig(
            'Cumulative Proportion of Variance Explained foe label_{}'.format(file_name))
        #plt.show()
        plt.close()
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
        fig.show()
        fig = px.scatter_3d(
            components, x=0, y=1, z=2, color=y_label,
            title=f'Total Explained Variance: {total_var:.2f}%',
            labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
        )
        fig.show()



    def RF_model(self,X_train, X_test, y_train, y_test,file_name="",i:int=0):
        """

        :param X_train:
        :param X_test:
        :param y_train:
        :param y_test:
        :return: randomforest model
        """

        hyperparameters = {'n_estimators': [100, 50],
                           'max_features': [0.3,0.5]
                           }
        #n_jobs number of cores used for it
        #scoring is the Strategy to evaluate the performance of the cross-validated model on the test set
        rf_cv = GridSearchCV(RandomForestClassifier(random_state=0),
                             hyperparameters, scoring='roc_auc_ovr_weighted',
                             cv=5,
                             verbose=True,
                             n_jobs=8)

        rf_cv.fit(X_train, y_train)
        print(rf_cv.best_params_)
        y_probs = rf_cv.best_estimator.predict_proba(X_test)
        # keep probabilities for the positive outcome only
        yhat = y_probs[:, 1]
        # calculate roc curves
        fpr, tpr, thresholds = roc_curve(y_test, yhat)
        # calculate the g-mean for each threshold
        gmeans = sqrt(tpr * (1 - fpr))
        # locate the index of the largest g-mean
        ix = argmax(gmeans)
        print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
        # plot the roc curve for the model
        plt.plot([0, 1], [0, 1], linestyle='--', label='No Skill')
        plt.plot(fpr, tpr, marker='.', label='Logistic')
        plt.scatter(fpr[ix], tpr[ix], marker='o', color='black',
                       label='Best')
        # axis labels
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend()
        # show the plot
        plt.savefig("ROC_curve{}".format(file_name))
        plt.show()
        plt.close()
        y_pred = rf_cv.best_estimator_.predict(X_test)
        cm = confusion_matrix(y_test, y_pred)
        cm_display = ConfusionMatrixDisplay(cm).plot()
        cm_display.figure_.savefig('cm_cv_best_parameters{}.png'.format(file_name), dpi=300)
        plt.close()
        #
        plt.figure()
        sns.lineplot(y=rf_cv.cv_results_["mean_test_score"],
                     x=rf_cv.cv_results_['param_max_features'].data,
                     hue=rf_cv.cv_results_['param_n_estimators']).set_xticks(
            range(22))
        plt.xlabel("max features")
        plt.ylabel("Accuracy (mean 10-fold CV)")
        plt.title(
            "Accuracy with different estimators and features for Amacrine RF model")
        plt.show()
        #use the best parameters from cross validation redo RF
        # rf = RandomForestClassifier(random_state=i,
        #                             n_estimators=100,max_features=0.5,
        #                             oob_score=True, n_jobs=6)
        # rf.fit(X_train, y_train)
        # y_pred = rf.predict(X_test)
        # print("Train/test accuracy: {}/{}".format(
        #     accuracy_score(rf.predict(X_train), y_train),
        #     accuracy_score(rf.predict(X_test), y_test)))
        #
        # print()
        # cm = confusion_matrix(y_test, y_pred)
        # cm_display = ConfusionMatrixDisplay(cm).plot()
        # cm_display.figure_.savefig('cm_{}.png'.format(file_name), dpi=300)
        # plt.title("RF confusion matrix with best parameters")
        # cm_display.figure_.savefig(
        #     'cm_cv_best_parameters{}.png'.format(file_name), dpi=300)
        # plt.show()
        #
        # fi = pd.DataFrame(data=rf.feature_importances_, index=X_train.columns,
        #                   columns=['Importance']) \
        #     .sort_values(by=['Importance'], ascending=False)
        #
        # And visualize
        sns.barplot(data=fi.head(20), x="Importance", y=(fi.head(20)).index).set_title(
            "feature importance for RF model")
        #fig = ax.get_figure()
        plt.savefig('feature_importance{}.png'.format(file_name))
        plt.show()
        # #ROC curve
        # fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        # self.show_roc(y_test, y_pred,  axes[1, 1],file_name)
        # plt.savefig()
        # show()

        # #save model
        # filename = 'data/model/rf_test_model_cv{}'.format(file_name)
        # joblib.dump(rf_cv, filename)
        # return rf_cv
    def show_roc(self,y_test,y_pred,ax,file_name):
        from sklearn.metrics import precision_score, accuracy_score, plot_roc_curve, RocCurveDisplay, \
            ConfusionMatrixDisplay, roc_auc_score, roc_curve
        fpr_rod, tpr_rod, _ = roc_curve(y_test, y_pred)
        auc_rod = roc_auc_score(y_test, y_pred)
        RocCurveDisplay(fpr=fpr_rod, tpr=tpr_rod, roc_auc=auc_rod,
                        estimator_name='Random forest prediction').plot(
            ax=ax)
        ax.set_title("Random forest on {}".format(file_name), fontsize=15)
        fig.suptitle(' ROC Curves at 0.5 threshold', fontsize=16)
        return fig,ax
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
        fingerprint_df=fingerprint_df_preapare(substrate_df,num_bits=num_bits,radius=3)
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

    def group_by_site(self,all_dataset=""):
        if all_dataset=="":
            with open('data/seq_smiles_all', 'rb') as file1:
                all_dataset = dill.load(file1)
        #add methyl_type
        for index in all_dataset.index:
            try:
                rec_sites = (all_dataset.loc[index,"reactant_site"]).split(",")
                sites=set()
                for i in rec_sites:
                    atom = list(i)
                    sites.add("".join(x for x in atom if x.isalpha()))
                print("".join(list(sites)))
                sites="_".join(list(sites))
                #only return alpha-->methyl atom symbol
                #["".join(x for x in rec_sites if x.isalpha())]
                all_dataset.loc[index, "methyl_type"] =sites
            except:
                all_dataset.loc[index, "methyl_type"] = 'NA'
        #save all data with methyl type
        with open('data/seq_smiles_all', 'wb') as dill_file:
            dill.dump(all_dataset,dill_file)
        seprate_dataset = all_dataset.groupby(by=["methyl_type"])
        for groups in seprate_dataset.groups:
            sub_df = seprate_dataset.get_group(groups)
            print(sub_df)
            sub_atom_object_dictionary = {}
            group = sub_df["methyl_type"].unique()
            print(group)
            for index in sub_df.index:
                sub_atom_object_dictionary[index]=sub_df.loc[index,"reactant_site"]
            print(sub_atom_object_dictionary)
            with open("data/group_data/seq_smiles_all{}".format(group), "wb") as dill_file:
                dill.dump(sub_df, dill_file)
            with open("data/group_data/diction_atom_all{}".format(group), "wb") as dill_file:
                dill.dump(sub_atom_object_dictionary, dill_file)



    def multidata_predict(self,dataset=None,inputmodel = None):
        if (dataset == None) or (inputmodel == None):
        #then do with test data
            loded_data=pd.read_excel("data/prediction_test.xlsx", header=0, index_col=0)
            print(loded_data)
            methyl_site_atom= predict(loded_data,num_bits=1024)
            print(methyl_site_atom)


def main():

    #run with command line
    # rh_file = argv[1]
    #run with pycharm




    data_with_site_drop_du = copy.deepcopy(data_with_site).drop_duplicates(['main_sub'])

    #
    # #read manual_data
    # parse_data.read_mannual_data()
    # with open("data/mannual_data", "rb") as dill_file:
    #     manual_data = dill.load(dill_file)
    # #save csv file to for check
    # with open("data/methyl_site_dictionary", "rb") as dill_file:
    #     methyl_site_dictionary = dill.load(dill_file)
    # #drop NA
    # indexNames = manual_data[manual_data['mainsub_mol'] == 'NA'].index
    # # Delete these row indexes from dataFrame
    # manual_data.drop(indexNames, inplace=True)
    # print(manual_data)

    # manual_fg = save_fingerprints_to_dataframe(manual_data,methyl_site_dictionary,2048,3,file_name="manual_2048")
    # print(manual_fg)
    # duplicate_1_class(manual_fg)


    #X = pd.read_csv("data/input_dataframe_withoutstructure_2048_drop_duplicate_all.csv",header=0,index_col=0)
    #X = duplicate_1_class(X,2)

    # merge manual data and rhea data

    #only use substrate then drop the duplicate


    #X_train, X_test, y_train, y_test = prepare_train_teat_data(X)

    #down sampling
    # from imblearn.under_sampling import TomekLinks
    # undersample = TomekLinks()
    # X_train, y_train = undersample.fit_resample(X_train, y_train)
    # #train RF model
    # model = RF_model(X_train, X_test, y_train, y_test,"_2048_drop_duplicate_double_true_label")
    # #multidata_predict()

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