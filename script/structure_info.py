#!/usr/bin/env python3
"""
need Numpy version equal or less than 1.21
"""
import copy
from urllib.request import urlopen

import dill
import pandas as pd
import prody as pro
from geometricus import MomentInvariants, SplitType
from geometricus import GeometricusEmbedding
import numpy as np
import umap
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import GroupShuffleSplit
from Model_class import Model_class

def download_pdb_structure_from_alphafola(input_file="",entry_pdb="../autodata/entry_pdb.xlsx"):

    input_df = pd.read_csv(input_file,header=0,index_col=0)
    entry_pdb_df=pd.read_excel(entry_pdb,header=0,index_col=0)
    input_df["structure"]=pd.DataFrame(
        len(input_df.index) * [0]).astype('object')
    input_df["invariants_kmer"]=pd.DataFrame(
        len(input_df.index) * [0]).astype('object')
    input_df["invariants_radius"]=pd.DataFrame(
        len(input_df.index) * [0]).astype('object')
    # # print(input_df.loc[2570,"structure"].dtype)
    # pdb_structure = pro.parsePDB(
    #     "../autodata/pdb_structure_from_alphafold/{}.pdb".format("A7HTX8"))
    # print(pdb_structure)
    #
    # input_df.loc[2570, "structure"] = pdb_structure
    # print(input_df)
    # print(entry_pdb_df)
    file1 = open(
        "../autodata/pdb_structure_from_alphafold/missing_structure.txt", "w")
    pdbs = []
    invariants_kmer = []
    invariants_radius = []
    missing_structure = []
    for index in input_df.index:
        print(index)
        entry = input_df.loc[index,"Entry"]
        #if there are structure in pdb database
        if entry in entry_pdb_df.index:
            try:

                pdb_id = entry_pdb_df.loc[entry,"PDB"].split(";")[0]
                print(pdb_id)
                pdb_structure = pro.parsePDB("../autodata/pdb_structure_from_pdb/{}.pdb".format(entry))
                pdbs.append(pdb_structure)
                invariants_kmer.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.KMER,
                                                          split_size=16))

                invariants_radius.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.RADIUS,
                                                          split_size=10))
                print(pdb_structure)
            except:
                pdb = (entry_pdb_df.loc[entry,"PDB"]).split(";")[0]
                url = "https://files.rcsb.org/download/{}.pdb".format(pdb)
                # Download from URL
                with urlopen(url) as webpage:
                    content = webpage.read()
                # Save to file
                with open("../autodata/pdb_structure_from_pdb/" + str(entry) + ".pdb",
                          'wb') as download:
                    download.write(content)
                print("save ../autodata/pdb_structure_from_pdb/" + str(entry) + ".pdb")
                # parase pdb and add to list
                pdb_structure = pro.parsePDB("../autodata/pdb_structure_from_pdb/{}.pdb".format(entry))
                pdbs.append(pdb_structure)
                invariants_kmer.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.KMER,
                                                          split_size=16))

                invariants_radius.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.RADIUS,
                                                          split_size=10))

        else:
            try:
                pdb_structure=pro.parsePDB(
                    "../autodata/pdb_structure_from_alphafold/{}.pdb".format(
                        entry))
                pdbs.append(pdb_structure)
                invariants_kmer.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.KMER,
                                                          split_size=16))

                invariants_radius.append(
                    MomentInvariants.from_prody_atomgroup(entry,
                                                          pdb_structure,
                                                          split_type=SplitType.RADIUS,
                                                          split_size=10))
            except:
                #otherwise download from alphafold
                try:
                    url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v3.pdb".format(
                        entry)
                    # Download from URL
                    with urlopen(url) as webpage:
                        content = webpage.read()
                    # Save to file
                    with open("../autodata/pdb_structure_from_alphafold/" + str(entry) + ".pdb",
                              'wb') as download:
                        download.write(content)
                    print("save ../autodata/pdb_structure_from_alphafold/" + str(entry) + ".pdb")
                    #parase pdb and add to list
                    pdb_structure = pro.parsePDB(
                        "../autodata/pdb_structure_from_alphafold/{}.pdb".format(
                            entry))
                    pdbs.append(pdb_structure)
                    invariants_kmer.append(
                        MomentInvariants.from_prody_atomgroup(entry,
                                                              pdb_structure,
                                                              split_type=SplitType.KMER,
                                                              split_size=16))

                    invariants_radius.append(
                        MomentInvariants.from_prody_atomgroup(entry,
                                                              pdb_structure,
                                                              split_type=SplitType.RADIUS,
                                                              split_size=10))
                except:
                    #save missing structure to file
                    file1.write("{}\n".format(entry))
                    print("{} do not have predicted structure in alphafold database".format(entry))
                    pdbs.append("NA")
                    invariants_kmer.append("NA")
                    invariants_radius.append("NA")
    file1.close()
    input_df["structure"] = pdbs


def geometrius_embedding(invariants_kmer,invariants_radius,entry):

        kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer,
                                                             resolution=2.,protein_keys=entry)
        radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius,
                                                               resolution=2.,protein_keys=entry)

        print(kmer_embedder )#2D
        print(radius_embedder)
        return kmer_embedder,radius_embedder

def merge_structure_embedding_to_input(input_df=""):

    #  with open("invariants_kmer", "rb") as dillfile1:
    #     invariants_kmer = dill.load(dillfile1)
    # dillfile1.close()
    #
    # with open("invariants_radius", "rb") as dillfile2:
    #     invariants_radius=dill.load(dillfile2)

    # dillfile2.close()
    # input_df["invariants_kmer"] = invariants_kmer
    # input_df["invariants_radius"] = invariants_radius
    with open(
            "../autodata/input_data/active_site/with_structure/PF08241_bit_score15_coverage0.7_ACS_bit128_3_remove_redundant_with_structure",
            "rb") as dillfile:
        input_df = dill.load(dillfile)
    print(input_df)

    # train_test split
    splitter = GroupShuffleSplit(test_size=0.8, n_splits=1, random_state=0)
    split = splitter.split(input_df, groups=input_df["molecular_id"])
    train_inds, test_inds = next(split)
    train = input_df.iloc[train_inds]
    test = input_df.iloc[test_inds]

    # X_train = train[list(train.columns)[:-2]]
    X_train = (copy.deepcopy(train)).drop(columns=["label"])
    Y_train = train["label"]
    # X_test = test[list(test.columns)[:-2]]
    X_test = (copy.deepcopy(test)).drop(columns=["label"])
    Y_test = test["label"]

    train_invariants_kmer = X_train["invariants_kmer"]
    train_invariants_radius = X_train["invariants_radius"]
    # embeddding for training data
    train_kmer_embedder, train_radius_embedder = geometrius_embedding(
        train_invariants_kmer, train_invariants_radius,entry=X_train["Entry"])
    shapemers_k_mer = train_kmer_embedder.shapemer_keys
    shapemers_radius_mer = train_radius_embedder.shapemer_keys
    # use the same shapmer for test embedding
    test_invariants_kmer = X_test["invariants_kmer"]
    test_invariants_radius = X_test["invariants_radius"]
    # test_kmer_embedder = GeometricusEmbedding.from_invariants(
    #     test_invariants_kmer,
    #     resolution=2., shapemer_keys=shapemers_k_mer,protein_keys=X_test["Entry"])
    # test_radius_embedder = GeometricusEmbedding.from_invariants(
    #     test_invariants_radius,
    #     resolution=2., shapemer_keys=shapemers_radius_mer,protein_keys=X_test["Entry"])
    test_kmer_embedder=train_kmer_embedder.embed(test_invariants_kmer, X_test["Entry"])
    print(len(train_kmer_embedder.embedding))
    print(train_kmer_embedder.embedding)
def main():
    merge_structure_embedding_to_input()
    #download_pdb_structure_from_alphafola(input_file="../autodata/input_data/active_site/PF08241_bit_score15_coverage0.7_ACS_bit128_3_remove_redundant.csv")

if __name__ == "__main__":
    main()