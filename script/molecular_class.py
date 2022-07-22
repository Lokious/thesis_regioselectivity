"""
class for saving the properties for each molecular
"""
from pikachu.general import draw_smiles
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display, Image
from rdkit.Chem import Draw
from rdkit import Chem
from PIL import Image
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
        else:
            i = index+1
            mol = Chem.MolFromSmiles(smile)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + i)
                # save the index in isotope just for keeeping the index for later use
                atom.SetIsotope(atom.GetIdx() + i)
                index = atom.GetIdx() + i
            return mol,index


class reaction ():
    def __init__(self,substrates="", products="",rxn_object=None):
        self.substrates = substrates
        self.products = products
        self.rxn_object = rxn_object
        self.mol_product = None
        self.mol_substrate = None

    def get_reaction_sites(self,rxn_object=""):
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
        with open("test.png", 'wb') as handle:
            handle.write(img)

    def perform_reaction(self):
        """
        perform reaction to get product's mol object and save, perform for each reaction

        :return:
        """
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
                Draw.ShowMol(mol_product, size=(600, 600))
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

        Draw.ShowMol(mol_substrate, size=(600, 600))
        Draw.ShowMol(mol_product, size=(600, 600))
        for atom in mol_substrate.GetAtoms():
            atom.SetAtomMapNum(atom.GetIsotope())
            #atom.SetIsotope(0)

        reactant_atoms = []
        pro_atom_dic= {}
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
        for index in reactant_atoms:
            atom_methyl = mol_product.GetAtomWithIdx(index)
            for bond in atom_methyl.GetBonds():
                atom_1, atom_2 = bond.GetBeginAtom(), bond.GetEndAtom()
                if atom_1 == atom_methyl:
                    atoms_list.append(atom_2)
                else:
                    atoms_list.append(atom_1)
        for atom in atoms_list:
            print(atom.GetIdx())
            print(atom.GetAtomMapNum())
            print(atom.GetIsotope())
        return atoms_list


def main():
    unittest.main()

if __name__ == "__main__":
    main()

