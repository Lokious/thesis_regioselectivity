"""
class for saving the properties for each molecular
"""
import pandas as pd
#from pikachu.general import draw_smiles
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display, Image
from rdkit.Chem import Draw,DataStructs
from rdkit import Chem
from PIL import Image
import numpy as np
import dill
import unittest

class molecular ():
    def __init__(self,chembi="",rhea_comp="",inchkey="",smiles=''):
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
    def calculate_smile(self):
        smile = Chem.MolToSmiles(mol)
        self.smiles = smile
    def mol_with_atom_index(self,smile="",mol_object = None, index ={}):
        """Return mol object with index, input could be simles or mol

        smile: string, smile of an molecular from substrates
        mol_object: create by rdkit or read from file
        index: dictionary, key is atom type<O,N,...> value is the largest index for this atom type plus 1
        """

        #if the input has mol_object
        if mol_object:
            #set i here inorder to make the atom from all substrates has different mapnumber
            atom_type = index
            for atom in mol_object.GetAtoms():

                atom_sy = atom.GetSymbol()
                if atom_sy not in atom_type.keys():
                    # makesure no repeate number
                    atom_type[atom_sy] = 0
                # print(atom_type[atom_sy])
                # print(atom.GetSymbol())
                atom.SetAtomMapNum(atom.GetIdx())
                # save the index in isotope just for keeeping the index for later use
                atom.SetIsotope(atom_type[atom_sy])
                atom_type[atom_sy] += 1

            return mol_object,atom_type
        #if no mol_object input, but have smile input
        elif smile:

            mol = Chem.MolFromSmiles(smile)
            if mol:
                atom_type = index
                for atom in mol.GetAtoms():

                    if atom.GetSymbol() not in atom_type.keys():
                        #makesure no repeate number for an atom type
                        atom_type[atom.GetSymbol()] = 0
                    atom.SetAtomMapNum(atom.GetIdx())
                    # save the index in isotope just for keeeping the index for later use
                    atom.SetIsotope(atom_type[atom.GetSymbol()])
                    atom_type[atom.GetSymbol()] += 1
                return mol,atom_type
            else:
                print("warning:cannot due with this smile:{}".format(smile))
                return None, index
        else:
            print("missing input")

    def create_fingerprint_mol(self, substrate_molecular: Chem.Mol, num_bits: int = 2048,
        radius: int = 3)->np.array:
        """
        This function is to create
        :param substrate_molecular:
        :param num_bits:
        :param radius:
        :return:
        """
        #sanitize molecular
        Chem.SanitizeMol(substrate_molecular,
                         sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        rdmolops.SanitizeFlags.SANITIZE_NONE
        #inisilize a numpy array for molecular fingerprint
        bit_fingerprint_mol = np.zeros((0,),
                                   dtype=int)  # (one dimention, 0 is number of rows)

        #returns an RDKit vector object.

        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(substrate_molecular, radius,
                                                                  num_bits)

        # convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_mol)

        return bit_fingerprint_mol

    def create_fingerprint_atom(self, substrate_molecular: Chem.Mol,atom_object: Chem.Atom, num_bits: int = 2048,
        radius: int = 3)->np.array:
        atom_index = atom_object.GetIdx()
        atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(substrate_molecular, radius,
                                                                 atom_index)
        atom_map = {}
        submol = Chem.PathToSubmol(substrate_molecular, atom_environment, atomMap=atom_map)
        #sanitize molecular
        Chem.SanitizeMol(submol,
                         sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        rdmolops.SanitizeFlags.SANITIZE_NONE
        #inisilize a numpy array for molecular fingerprint
        bit_fingerprint_atom = np.zeros((0,),
                                   dtype=int)  # (one dimention, 0 is number of rows)

        #returns an RDKit vector object.

        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(submol, radius,
                                                                  num_bits)

        # We convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_atom)

        return bit_fingerprint_atom

class reaction ():
    def __init__(self,substrates="", products="",rxn_object=None):
        self.substrates = substrates
        self.products = products
        self.rxn_object = rxn_object
        self.mol_product = None
        self.mol_substrate = None

    def get_reaction_sites(self,rxn_object="",file_name= None ):
        """
        ###not use anymore, just to save the picture of reaction####

        :param rxn_object:
        :return:
        """
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

    def perform_reaction(self):
        """
        perform reaction to get product's mol object and save, perform for each reaction
        ###This function is not use anymore, it is before i want to get the atom map between substrate and product
        but you should set the some atom with same mapnumber to let the program know they are the same atom###

        :return:
        """
        # Draw.ShowMol(self.substrates, size=(600, 600))
        # Draw.ShowMol(self.products, size=(600, 600))
        #link the main substrate molecular and main product, to build a reaction template
        react = AllChem.ReactionFromSmarts(
            self.substrates + ">>" + self.products, useSmiles=True)
        mol_substrate = Chem.MolFromSmiles(self.substrates)

        self.mol_substrate = mol_substrate
        for atom in mol_substrate.GetAtoms():
            atom.SetIsotope(0)
        # print(self.products)
        # print(self.substrates)
        react.Initialize()
        print(self.products)
        print(self.substrates)
        # There should be only one product molecular,
        # the RunReactants function just for keep the atom map to find the reactant atom
        for result in react.RunReactants((mol_substrate,)):
            # perform reaction and get the product
            for mol_product in result:
                #Draw.ShowMol(mol_product, size=(600, 600))
                for atom in mol_product.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIsotope())
                    print(atom.GetAtomMapNum())
                    # atom.SetIsotope(0)
                result1 = Chem.MolToSmiles(mol_product)
        self.mol_product = mol_product

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

    def fingerprint_similiarity(self,mol1_fprint,mol2_fprint,mol1:Chem.Mol,mol2:Chem.Mol):

        # First we check if the dimensions of both fingerprints are correct.
        if len(mol1_fprint.shape) != 1:
            raise ValueError(f"expected dimensionality (N,) for `first_fingerprint`, got: {mol1_fprint.shape}")
        if len(mol2_fprint.shape) != 1:
            raise ValueError(f"expected dimensionality (N,) for `second_fingerprint`, got: {second_fingerprint.shape}")

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


    def main_substrate(self,subs,pros):
        """
        This is the function for find the substrate which is methylated among all and the product

        :param subs: list of substrate smiles from one reaction
        :param pros: list of product smiles from the same reaction
        :return:the substrate which is methylated among all and the corresponding product
        """
        sim_dictionary = {}
        mol_object = molecular()
        reaction_object = reaction()
        for i,mol1_smile in enumerate(subs):
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

                sim_dictionary[
                    (i, j)] = reaction_object.fingerprint_similiarity(sub_fg,pro_fg,mol1=mol1,
                                                                      mol2=mol2)
        similarity_list_top_2 = list(
            sorted(sim_dictionary.items(),reverse=True,key=lambda item: item[1]))[:2]
        for key in similarity_list_top_2:
            i = key[0][0]
            j= key[0][1]

            #this function assumed the similarity of fingerprint betwween molecular before methylation
            #and after methylation should be higher than this molecular with other molecular
            #the methyl donor molecular will become smaller after reaction
            if (len(pros[j])>len(subs[i])) and (pros[j] != 0) and (subs[i] != 0 ):
                return [subs[i],pros[j]]
            else:
                continue
class Testreaction_class(unittest.TestCase):
    def test0_get_reactant_atom(self):
        reaction_obj =reaction(substrates="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1OH:7])[4cH:4][2cH:2]1",
                               products="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1O:7][11CH3:13])[4cH:4][2cH:2]1")
        atom_list, index_list, mainsub_mol = reaction_obj.get_reactant_atom()
        self.assertEqual(index_list[0],"O7:1")
def main():
    unittest.main()

if __name__ == "__main__":
    main()

