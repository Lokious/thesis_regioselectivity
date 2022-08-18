"""
class for saving the properties for each molecular  # this docstring is not clear to me
"""
import copy
import typing
import pandas as pd
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display, Image
from rdkit.Chem import Draw,DataStructs
from rdkit import Chem
from PIL import Image
import numpy as np
import dill
import unittest


class Molecule:  # classes are written in CamelCase
    def __init__(self, chembi="", rhea_comp="", inchkey="", smiles=""):
        self.chembi = chembi
        self.rhea_comp = rhea_comp
        self.inchkey = inchkey
        self.smiles = smiles

    def get_chembi(self):
        return self.chembi
    
    def get_rhea_comp(self):
        return self.rhea_comp
    
    def get_inchkey(self):
        return self.inchkey
    
    def get_smiles(self):
        return self.smiles

    def mol_with_atom_index(self, input: typing.Optional[typing.Union[str, Chem.Mol]] = None, index: dict = dict()):
        """Return mol object with index, input could be smiles or mol #  it is not clear to me from the docstring what the returned index is

        smile: string, smile of an molecular from substrates
        mol_object: create by rdkit or read from file
        index: dictionary, key is atom type<O,N,...> value is the largest index for this atom type plus 1
        """
        atom_type = index
        if isinstance(input, str): 
            # evaluate if SMILES string is valid...
            mol = Chem.MolFromSmiles(input)
            if mol:

                for atom in mol.GetAtoms():

                    if atom.GetSymbol() not in atom_type.keys():
                        # make sure no repeat number for an atom type  # repeat number?
                        atom_type[atom.GetSymbol()] = 0
                    atom.SetAtomMapNum(atom.GetIdx())
                    # save the index in isotope just for keeeping the index for later use
                    atom.SetIsotope(atom_type[atom.GetSymbol()])
                    atom_type[atom.GetSymbol()] += 1

                return mol, atom_type
            else:
                print(f"warning: can not deal with this SMILES: {input}")  # what does this message mean exactly? I know this is draft code, but make sure the English is correct
                return None, index
        elif isinstance(input, Chem.Mol):
            # do something with RDKit Mol object...
                # set i here in order to make the atom from all substrates has different mapnumber
            for atom in input.GetAtoms():
                atom_symbol = atom.GetSymbol()  # use as little abbreviations as possible, makes code easier to understand
                if atom_symbol not in atom_type.keys():
                    # make sure no repeat number  # watch the spelling
                    atom_type[atom_symbol] = 0
                atom.SetAtomMapNum(atom.GetIdx())
                # save the index in isotope just for keeping the index for later use
                atom.SetIsotope(atom_type[atom_symbol])
                atom_type[atom_symbol] += 1

            return input, atom_type
        else:
            print("missing input")

    def create_fingerprint_mol(self, substrate_molecular: Chem.Mol, num_bits: int = 2048,
        radius: int = 3) -> np.array:
        """
        This function is to create  # to create what? :)
        
        :param substrate_molecular:
        :param num_bits:
        :param radius:
        :return:
        """
        #sanitize molecular
        Chem.SanitizeMol(
            substrate_molecular,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
        rdmolops.SanitizeFlags.SANITIZE_NONE
        # initialize a numpy array for molecular fingerprint
        bit_fingerprint_mol = np.zeros(
            (0,),
            dtype=int)  # (one dimention, 0 is number of rows)

        #returns an RDKit vector object.

        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(
            substrate_molecular, 
            radius,
            num_bits
        )

        # convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_mol)

        return bit_fingerprint_mol

    def create_fingerprint_atom(
        self, 
        substrate_molecular: Chem.Mol,
        atom_object: Chem.Atom, 
        num_bits: int = 2048,
        radius: int = 3
    ) -> np.array:
        """
        documentation?
        """
        atom_index = atom_object.GetIdx()
        atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
            substrate_molecular, 
            radius,
            atom_index
        )
        atom_map = {}
        submol = Chem.PathToSubmol(substrate_molecular, atom_environment, atomMap=atom_map)
        # sanitize molecule
        Chem.SanitizeMol(submol,
                         sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        rdmolops.SanitizeFlags.SANITIZE_NONE
        #inisilize a numpy array for molecular fingerprint
        bit_fingerprint_atom = np.zeros(
            (0,),
            dtype=int  # (one dimention, 0 is number of rows)
        )

        #returns an RDKit vector object.
        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(
            submol, 
            radius,
            num_bits
        )
        
        # We convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_atom)

        return bit_fingerprint_atom

    
class Reaction:
    def __init__(self, substrates="", products="", rxn_object=None):
        self.substrates = substrates
        self.products = products
        self.rxn_object = rxn_object
        self.mol_product = None
        self.mol_substrate = None

    def perform_reaction(self,rxn_object="",file_name= None ):
        """
        ###not use anymore, just to save the picture of reaction####

        :param rxn_object:
        :return:
        """
        raise RuntimeError("function `get_reaction_sites()` is deprecated")
        
        rxn_object.Initialize()
        
        try:
            reacting_atom_set = rxn_object.GetReactingAtoms()
            print(reacting_atom_set)
            img = Draw.ReactionToImage(rxn_object, returnPNG=True, subImgSize=(600,600))
            img.drawOptions()
            with open("{}.png".format(file_name), 'wb') as handle:
                handle.write(img)
            handle.close()
        except:
            print("canot run reaction{}".format(file_name))

    def get_reaction_sites(self,product_smiles:str="",substrate_smiles:str=""):
        """
        finding methylation atoms by removing -CH3 group and compare the similarity between substrate and product
        :param product_smiles:sting,the SMILES of Molecule after methylation
            substrate_smiles:sting,the SMILES of Molecule before methylation
        :return:
        """
        #remove AtomMapNum and Isotope then they will not affect similarity
        product_mol = Chem.MolFromSmiles(product_smiles)
        for atom in product_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)
        substrate_mol = Chem.MolFromSmiles(substrate_smiles)
        for atom in substrate_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)
        #Draw.ShowMol(product_mol,(600,600))

        #get editable mol for removing methylation group
        mol = Chem.EditableMol(product_mol)
        atom_index_list = []
        similiarity = 0
        possiable_list=[]
        for atom in product_mol.GetAtoms():
            #do not consider C-methylation here
            if atom.GetSymbol() !="C":
                neighbours = atom.GetNeighbors()
                for neighbour_atom in neighbours:
                    if ((neighbour_atom.GetSymbol()=="C") and (neighbour_atom.GetDegree()==1)):
                        #print(atom.GetIdx(),atom.GetSymbol())
                        #mol_save_status = mol
                        try:
                            #save possiable atoms to a list, if remove one do not get 100% similarity substrate
                            #the list will be use later
                            mol = Chem.EditableMol(product_mol)
                            mol.RemoveAtom(neighbour_atom.GetIdx())
                            possiable_list.append(neighbour_atom.GetIdx())
                            # Draw.ShowMol(mol.GetMol(), (600, 600))
                            # Draw.ShowMol(substrate_mol, (600, 600))
                        except:
                            continue
                        similiarity = DataStructs.FingerprintSimilarity(
                            Chem.RDKFingerprint(mol.GetMol()),
                            Chem.RDKFingerprint(substrate_mol))
                        print(similiarity)
                        if similiarity==1:
                            atom_index_list.append(
                                atom.GetSymbol() + ":"+str(
                                    atom.GetIdx()))
                            mol_remove_methylation = mol.GetMol()
                            print("check")
                            return product_mol,Chem.MolToSmiles(mol_remove_methylation),atom_index_list,"Pass_Check"
                        else:
                            #mol=mol_save_status
                            continue
        #remove multiple methylation group
        from itertools import combinations
        index_combination_list=list(combinations(possiable_list, 2))
        index_combination_list += list(combinations(possiable_list, 3))
        print(list(combinations(possiable_list, 2)))
        for items in index_combination_list:
            mol=Chem.EditableMol(product_mol)
            remove_atom=sorted(list(items),reverse=True)
            #remove from larger index
            for i in remove_atom:
                mol.RemoveAtom(i)
            #Draw.ShowMol(mol.GetMol(), (600, 600))
            similiarity = DataStructs.FingerprintSimilarity(
                Chem.RDKFingerprint(mol.GetMol()),
                Chem.RDKFingerprint(substrate_mol))
            print(similiarity)
            if similiarity ==1:
                atom = (product_mol.GetAtomWithIdx(i)).GetNeighbors()[0]
                atom_index_list.append(
                    atom.GetSymbol() + ":"+str(
                        atom.GetIdx()))
                print(atom.GetSymbol() + ":"+str(
                        atom.GetIdx()))
                mol_remove_methylation = mol.GetMol()
                #Draw.ShowMol(mol.GetMol(), (600, 600))
                print("check")
                return product_mol, Chem.MolToSmiles(
                    mol_remove_methylation), atom_index_list, "Pass_Check"

        mol_remove_methylation=mol.GetMol()

        #print(similiarity)
        #raise RuntimeError("function `get_reaction_sites()` is deprecated")

        return product_mol,Chem.MolToSmiles(
            mol_remove_methylation), atom_index_list, "unCheck"
        
    def get_reactant_atom(self):
        """
        This function only consider the large substrate and product, it assume
        the smaller molecular as Methyl donor, H2O or coenzyme etc.

        :return: list of atom objects which is the regioselectivity site
        """
        #list of atom object
        atoms_list=[]
        atom_index_list = []
        #the main substrate and product to mol object, here the substrates and products only contains the main one
        mol_substrate = Chem.MolFromSmiles(self.substrates)
        mol_product = Chem.MolFromSmiles(self.products)

        sub_mapnumber_dict = {}
        pro_mapnumber_dict = {}

        #save the symbol and isotope as the key, and atom index as value in dictionary
        for atom in mol_substrate.GetAtoms():
            #isotope is what is set before
            sym_and_num = atom.GetSymbol()+str(atom.GetIsotope())
            index = atom.GetIdx()
            sub_mapnumber_dict[sym_and_num]=index
        for atom in mol_product.GetAtoms():
            sym_and_num = atom.GetSymbol()+str(atom.GetIsotope())
            index = atom.GetIdx()
            pro_mapnumber_dict[sym_and_num]=index

        for map_index in sub_mapnumber_dict.keys():
            atom_sub = mol_substrate.GetAtomWithIdx(sub_mapnumber_dict[map_index])
            if map_index in pro_mapnumber_dict.keys():
                atom_pro = mol_product.GetAtomWithIdx(pro_mapnumber_dict[map_index])
            else:
                continue
            #if atom neighbor increase and this atom is not carbon nor R group, then add it to methylation site list
            #now the problem here is the index is different in product and substrate sometime, then it will get wrong atom
            if len(atom_sub.GetNeighbors()) < len(atom_pro.GetNeighbors()):
                if (atom_sub.GetSymbol()!="C") and (atom_sub.GetSymbol()!="*") :
                    atoms_list.append(atom_sub)

        for atom in atoms_list:
            #str with symbol index and mapnumber(mapnumber is the same as substrate's mapnumber)
            atom_index_list.append((atom.GetSymbol() + str(atom.GetAtomMapNum())+":"+str(atom.GetIsotope())))
            # print(atom.GetAtomMapNum())
            # print(atom.GetIsotope())
        #print(atom_index_list)

        return atoms_list, atom_index_list,mol_substrate

    def fingerprint_similiarity(self,mol1_fprint="",mol2_fprint="",mol1:Chem.Mol=None,mol2:Chem.Mol=None):



        # First we check if the dimensions of both fingerprints are correct.
        if len(mol1_fprint.shape) != 1:
            raise ValueError(f"expected dimensionality (N,) for `first_fingerprint`, got: {mol1_fprint.shape}")
        if len(mol2_fprint.shape) != 1:
            raise ValueError(f"expected dimensionality (N,) for `second_fingerprint`, got: {mol2_fprint.shape}")

        # We also check if the lengths of both fingerprints are equal.
        if mol1_fprint.shape[0] != mol2_fprint.shape[0]:
            raise ValueError(
                f"first_fingerprint (num_bits: {mol1_fprint.shape[0]}) and "
                f"second_fingerrint (num_bits: {mol2_fprint.shape[0]}) do not "
                "have same length!"
            )

        tanimoto_similarity = (
            np.logical_and(mol1_fprint, mol2_fprint).sum() / (
                float(np.logical_or(mol1_fprint, mol2_fprint).sum())
            )
        )
        if mol1 and mol2:
            similarity = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2))
            return similarity
        return tanimoto_similarity

def main_substrate(subs, pros):
    """
    This is the function for find the substrate which is methylated among all and the product

    :param subs: list of substrate smiles from one reaction
    :param pros: list of product smiles from the same reaction
    :return:the substrate which is methylated among all and the corresponding product
    """
    sim_dictionary = {}
    mol_object = Molecule()
    reaction_object = Reaction()

    for i, mol1_smile in enumerate(subs):
        try:
            mol1 = Chem.MolFromSmiles(mol1_smile)
            sub_fg = mol_object.create_fingerprint_mol(substrate_molecular=mol1)
        except:
            return None

        for j,mol2_smile in enumerate(pros):
            try:
                mol2 = Chem.MolFromSmiles(mol2_smile)
                pro_fg = mol_object.create_fingerprint_mol(substrate_molecular=mol2)
            except:
                return None

            sim_dictionary[(i, j)] = reaction_object.fingerprint_similiarity(
                sub_fg,
                pro_fg,
                mol1=mol1,
                mol2=mol2
            )

    similarity_list_top_2 = list(
        sorted(
            sim_dictionary.items(),
            reverse=True,
            key=lambda item: item[1]
        )
    )[:2]

    for key in similarity_list_top_2:
        i = key[0][0]
        j= key[0][1]

        # the documentation below is not quite clear to me

        #this function assumed the similarity of fingerprint betwween molecular before methylation
        #and after methylation should be higher than this molecular with other molecular
        #the methyl donor molecular will become smaller after reaction
        if (len(pros[j]) > len(subs[i])) and (pros[j] != 0) and (subs[i] != 0 ):
            return [subs[i], pros[j]]
        else:
            continue
                
                
# class Testreaction_class(unittest.TestCase):
#     def test0_get_reactant_atom(self):
#         """
#         don't forget documentation for the test :)
#         """
#         reaction_obj =Reaction(
#             substrates="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1OH:7])[4cH:4][2cH:2]1",
#             products="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1O:7][11CH3:13])[4cH:4][2cH:2]1"
#         )
#         atom_list, index_list, mainsub_mol = reaction_obj.get_reactant_atom()
#         self.assertEqual(index_list[0], "O7:1")
#
#
def main():
    #unittest.main()
    reaction = Reaction()
    df = pd.read_csv("../data/seq_smiles_all_test_codefor_finding_methylsite.csv",header=0,index_col=0)
    df.dropna(inplace=True)
    df_save=copy.deepcopy(df)
    df_save["remove_methylation"]=pd.DataFrame(
            len(df_save.index) * [0]).astype('object')
    df_save["list_methylsite"]=pd.DataFrame(
            len(df_save.index) * [0]).astype('object')
    df_save["check"]=pd.DataFrame(
            len(df_save.index) * [0]).astype('object')
    for index in df.index:
        print(index)
        main_sub_smile = df.loc[index,"main_sub"]
        main_pro_smile = df.loc[index,"main_pro"]
        pro_mol,remove_methyl_smile,list_methylsite,check=reaction.get_reaction_sites(main_pro_smile,main_sub_smile)
        if check=="check":
            Draw.MolToFile(Chem.MolFromSmiles(remove_methyl_smile), "../pro_fig/{}_remove.png".format(index))
        Draw.MolToFile(pro_mol,
                       "../pro_fig/{}_pro.png".format(index))
        df_save.loc[index,"remove_methylation"]=remove_methyl_smile
        df_save.loc[index,"list_methylsite"]=",".join(list_methylsite)
        df_save.loc[index,"check"]=check
    print(df_save)
    df_save.to_csv("test_finding_methylation_site.csv")
if __name__ == "__main__":
    main()
