#!/usr/bin/env python3
"""
This works on server, seems due to environment problem of Biopython and prody on my pycharm
"""
from urllib.request import urlopen
import pandas as pd
import prody as pro
def download_pdb_structure_from_alphafola(input_file="",entry_pdb="../autodata/entry_pdb.xlsx"):

    input_df = pd.read_csv(input_file,header=0,index_col=0)
    entry_pdb_df=pd.read_excel(entry_pdb,header=0,index_col=0)
    input_df["structure"]=pd.DataFrame(
        len(input_df.index) * [0]).astype('object')

    print(entry_pdb_df)
    file1 = open(
        "../autodata/pdb_structure_from_alphafold/missing_structure.txt", "w")

    for index in input_df.index:
        entry = input_df.loc[index,"Entry"]
        #if there are structure in pdb database
        if entry in entry_pdb_df.index:
            try:
                pdbs=[]
                pdb_id=entry_pdb_df.loc[entry,"PDB"]
                pro.parsePDB("../autodata/pdb_structure_from_[db]/{}.pdb".format(pdb_id))
                input_df.loc[index,"structure"]=pdbs
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
        else:
            try:
                pdbs=[]
                pdb_id=entry_pdb_df.loc[entry,"PDB"].split(";")[0]
                pro.parsePDB(
                    "../autodata/pdb_structure_from_alphafold/{}.pdb".format(
                        pdb_id))
                input_df.loc[index,"structure"]=pdbs
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
                except:
                    #save missing structure to file
                    file1.write("{}\n".format(entry))
                    print("{} do not have predicted structure in alphafold database".format(entry))

def main():
    download_pdb_structure_from_alphafola(input_file="../autodata/input_data/input128fg_bond3_S_k_mer.csv")

if __name__ == "__main__":
    main()