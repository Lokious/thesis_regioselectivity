#!/usr/bin/env python3
"""
This works on server, seems due to environment problem of Biopython and prody on my pycharm
"""
from urllib.request import urlopen
import pandas as pd
import prody as pro
from geometricus import MomentInvariants, SplitType
from geometricus import GeometricusEmbedding
import numpy as np
import umap
import matplotlib.pyplot as plt
def download_pdb_structure_from_alphafola(input_file="",entry_pdb="../autodata/entry_pdb.xlsx"):

    input_df = pd.read_csv(input_file,header=0,index_col=0)
    entry_pdb_df=pd.read_excel(entry_pdb,header=0,index_col=0)
    input_df["structure"]=pd.DataFrame(
        len(input_df.index) * [0]).astype('object')
    # print(input_df.loc[2570,"structure"].dtype)
    pdb_structure = pro.parsePDB(
        "../autodata/pdb_structure_from_alphafold/{}.pdb".format("A7HTX8"))
    print(pdb_structure)
    invariants_kmer = [MomentInvariants.from_prody_atomgroup("A7HTX8", pdb_structure,split_type=SplitType.KMER,
                                                                split_size=16)]

    invariants_radius = [MomentInvariants.from_prody_atomgroup("A7HTX8", pdb_structure,
                                              split_type=SplitType.RADIUS,
                                              split_size=10)]
    kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer,
                                                         resolution=2.)
    radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius,
                                                           resolution=2.)

    print(kmer_embedder )#2D

    input_df.loc[2570, "structure"] = pdb_structure
    print(input_df)
    print(entry_pdb_df)
    file1 = open(
        "../autodata/pdb_structure_from_alphafold/missing_structure.txt", "w")

    for index in input_df.index:
        print(index)
        entry = input_df.loc[index,"Entry"]
        #if there are structure in pdb database
        if entry in entry_pdb_df.index:
            try:
                pdbs=[]
                pdb_id = entry_pdb_df.loc[entry,"PDB"].split(";")[0]
                print(pdb_id)
                pdb_structure = pro.parsePDB("../autodata/pdb_structure_from_pdb/{}.pdb".format(entry))
                print(pdb_structure)
                input_df.loc[index,"structure"] = pdb_structure
                print(input_df)
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
            print(input_df)
        else:
            print(input_df)
            try:
                pdbs=[]
                pdb_structure=pro.parsePDB(
                    "../autodata/pdb_structure_from_alphafold/{}.pdb".format(
                        entry))
                input_df.loc[index,"structure"]=pdb_structure
                print(input_df)
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
                except:
                    #save missing structure to file
                    file1.write("{}\n".format(entry))
                    print("{} do not have predicted structure in alphafold database".format(entry))
    print(input_df)
def main():

    download_pdb_structure_from_alphafola(input_file="../autodata/input_data/active_site/PF08241_bit_score15_coverage0.7_ACS_bit128_3_remove_redundant.csv")

if __name__ == "__main__":
    main()