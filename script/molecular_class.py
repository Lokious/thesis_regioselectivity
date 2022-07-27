"""
class for saving the properties for each molecular
"""
import pandas as pd
from pikachu.general import draw_smiles
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display, Image
from rdkit.Chem import Draw,DataStructs
from rdkit import Chem
from PIL import Image
import numpy as np
import dill
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
    def mol_with_atom_index(self,smile="",mol_object = None, index = 0):
        """Return mol object with index, input could be simles or mol"""
        if mol_object:
            #set i here inorder to make the atom from all substrates has different mapnumber
            i = index+1
            for atom in mol_object.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + i)
                # save the index in isotope just for keeeping the index for later use
                atom.SetIsotope(atom.GetIdx() + i)
                index = atom.GetIdx() + i
            return mol_object,index
        elif smile:
            i = index+1
            mol = Chem.MolFromSmiles(smile)
            if mol:
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIdx() + i)
                    # save the index in isotope just for keeeping the index for later use
                    atom.SetIsotope(atom.GetIdx() + i)
                    index = atom.GetIdx() + i
                return mol,index
            else:
                print("warning:cannot due with this smile:{}".format(smile))
                return None, index
        else:
            print("missing input")
    def create_fingerprint_mol(self, substrate_molecular: Chem.Mol, num_bits: int = 2048,
        radius: int = 3)->np.array:
        """

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

        # We convert the RDKit vetor object to a numpy array.
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
        Function for finding all reactanct atoms, which include all substrates
        and products in the reaction

        :param rxn_object:
        :return:
        """
        rxn_object.Initialize()
        reacting_atom_set = rxn_object.GetReactingAtoms()
        # print(reacting_atom_set)
        img = Draw.ReactionToImage(rxn_object, returnPNG=True, subImgSize=(600,600))
        img.drawOptions()
        with open("{}.png".format(file_name), 'wb') as handle:
            handle.write(img)
        handle.close()


    def perform_reaction(self):
        """
        perform reaction to get product's mol object and save, perform for each reaction

        :return:
        """
        # Draw.ShowMol(self.substrates, size=(600, 600))
        # Draw.ShowMol(self.products, size=(600, 600))
        print(self.substrates)
        print(self.products)
        react = AllChem.ReactionFromSmarts(
            self.substrates + ">>" + self.products, useSmiles=True)
        mol_substrate = Chem.MolFromSmiles(self.substrates)
        self.mol_substrate = mol_substrate
        # print(self.products)
        # print(self.substrates)
        react.Initialize()
        # There should be only one product molecular,
        # the RunReactants function just for keep the atom map to find the reactant atom
        for result in react.RunReactants((mol_substrate,)):
            # perform reaction and get the product
            for mol_product in result:
                #Draw.ShowMol(mol_product, size=(600, 600))
                for atom in mol_product.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIsotope())
                    # atom.SetIsotope(0)
                result1 = Chem.MolToSmiles(mol_product)
        self.mol_product = mol_product

    def get_reactant_atom(self):
        """
        This function only consider the large substrate and product, it assume
        the smaller molecular as Methyl donor, H2O or coenzyme etc.

        :return: list of atom objects which is the regioselectivity site
        """
        #perform the reaction to
        reaction.perform_reaction(self)

        mol_product = self.mol_product
        mol_substrate = self.mol_substrate
        # Draw.ShowMol(mol_substrate, size=(600, 600))
        # Draw.ShowMol(mol_product, size=(600, 600))
        for atom in mol_substrate.GetAtoms():
            atom.SetAtomMapNum(atom.GetIsotope())
            #atom.SetIsotope(0)

        reactant_atoms = []
        pro_atom_dic = {}
        sub_atom_dic = {}
        for atom in mol_product.GetAtoms():
            pro_atom_dic[atom] = atom.GetAtomMapNum()
        for atom in mol_substrate.GetAtoms():
            sub_atom_dic[atom] = atom.GetAtomMapNum()
        #list of tuple with (atom_object,AtomMapNum)
        pro_atom_list = list(sorted(pro_atom_dic.items(), key=lambda item: item[1]))
        sub_atom_list = list(
            sorted(sub_atom_dic.items(), key=lambda item: item[1]))

        for i in range(len(pro_atom_list)):

            atom_index_pro = pro_atom_list[i][0].GetIdx()
            radius = 1
            try:
                atom_index_sub = sub_atom_list[i][0].GetIdx()
                sub_atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(mol_substrate,
                                                                         radius,
                                                                         atom_index_sub)
            except:
                reactant_atoms.append(atom_index_pro)

            pro_atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
                mol_product,
                radius,
                atom_index_pro)
            atom_map = {}
            sub_mol = Chem.PathToSubmol(mol_substrate, sub_atom_environment,
                              atomMap=atom_map)
            pro_mol = Chem.PathToSubmol(mol_product, pro_atom_environment,
                              atomMap=atom_map)

            # Draw.ShowMol(sub_mol, size=(600, 600))
            # Draw.ShowMol(pro_mol, size=(600, 600))

        atoms_list=[]
        atom_index_list = []

        for index in reactant_atoms:
            atom_methyl = mol_product.GetAtomWithIdx(index)
            for bond in atom_methyl.GetBonds():
                atom_1, atom_2 = bond.GetBeginAtom(), bond.GetEndAtom()
                if atom_1 == atom_methyl:
                    atoms_list.append(atom_2)

                else:
                    atoms_list.append(atom_1)

        for atom in atoms_list:
            #str with symbol index and mapnumber(mapnumber is the same as substrate's mapnumber)
            atom_index_list.append((atom.GetSymbol() + str(atom.GetIdx())+":"+str(atom.GetAtomMapNum())))
            # print(atom.GetAtomMapNum())
            # print(atom.GetIsotope())
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
                #sim_dictionary[(i,j)]=reaction_object.fingerprint_similiarity(sub_fg,pro_fg)
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

def main():
    unittest.main()

if __name__ == "__main__":
    main()

