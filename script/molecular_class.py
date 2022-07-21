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

    def get_reaction_sites(self,rxn_object=""):
        """
        Function for finding all reactanct atoms, which include all substrates
        and products in the reaction

        :param rxn_object:
        :return:
        """
        rxn_object.Initialize()
        reacting_atom_set = rxn_object.GetReactingAtoms()
        print(reacting_atom_set)
        img = Draw.ReactionToImage(rxn_object, returnPNG=True)
        with open("test.png", 'wb') as handle:
            handle.write(img)

    def get_reactant_atom(self):
        """
        This function only consider the large substrate and product, it assume
        the smaller molecular as Methyl donor, H2O or coenzyme etc.

        :return:
        """
        react = AllChem.ReactionFromSmarts(
            self.substrates + ">>" + self.products, useSmiles=True)
        mol_substrate = Chem.MolFromSmiles(self.substrates)
        print(self.products)
        react.Initialize()
        #There should be only one product molecular,
        #the RunReactants function just for keep the atom map to find the reactant atom
        for result in react.RunReactants((mol_substrate,)):
            for mol_product in result:
                #Draw.ShowMol(mol2, size=(600, 600))
                for atom in mol_product.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIsotope())
                    #atom.SetIsotope(0)
                result = Chem.MolToSmiles(mol_product)
                print(result)
                print(result == self.products)
        for atom in mol_substrate.GetAtoms():
            atom.SetAtomMapNum(atom.GetIsotope())
            #atom.SetIsotope(0)
        for atom in mol_product.GetAtoms():
            atom_index = atom.GetIdx()
            radius = 1
            reactant_atoms = []
            try:
                sub_atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(mol_substrate,
                                                                         radius,
                                                                         atom_index)
            except:
                reactant_atoms.append(atom_index)
            pro_atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
                mol_product,
                radius,
                atom_index)

            atom_map = {}
            sub_submol_smiles = Chem.MolToSmiles(Chem.PathToSubmol(mol_substrate, sub_atom_environment, atomMap=atom_map))
            pro_submol_smiles = Chem.MolToSmiles(Chem.PathToSubmol(mol_product, pro_atom_environment,
                                       atomMap=atom_map))
            print(sub_submol_smiles)
            print(pro_submol_smiles)
            if sub_submol_smiles != pro_submol_smiles:
                reactant_atoms.append(atom_index)
                print(reactant_atoms)
        im = Draw.FingerprintEnv(mol_substrate,highlightAtomLists=reactant_atoms)
        im.show()
        #Draw.ShowMol(mol1, size=(600, 600))



def main():
    unittest.main()

if __name__ == "__main__":
    main()

