"""
class for saving the properties for each molecular
"""
from pikachu.general import draw_smiles
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
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
    def mol_with_atom_index(self,smile="",mol_object = None):
        """Return mol object with index, input could be simles or mol"""
        if mol_object:
            for atom in mol_object.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
                # save the index in isotope just for keeeping the index for later use
                atom.SetIsotope(atom.GetIdx() + 1)
            return mol_object
        else:
            mol = Chem.MolFromSmiles(smile)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
                # save the index in isotope just for keeeping the index for later use
                atom.SetIsotope(atom.GetIdx() + 1)
            return mol


class reaction ():
    def __init__(self,substrates="", products="",rxn_object=None):
        self.substrates = substrates
        self.products = products
        self.rxn_object = rxn_object

    def get_reaction_sites(self,rxn_object=""):
        rxn_object.Initialize()
        reacting_atom_set = rxn_object.GetReactingAtoms()
        print(reacting_atom_set)
        img = Draw.ReactionToImage(rxn_object, returnPNG=True)
        with open("test.png", 'wb') as handle:
            handle.write(img)

    def get_reactant_atom(self):

        react = AllChem.ReactionFromSmarts(
            self.substrates + ">>" + self.products, useSmiles=True)
        results = react.RunReactants((substrates,))[0]
        print(results == self.products)


def main():
    unittest.main()

if __name__ == "__main__":
    main()

