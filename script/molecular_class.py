"""
class for saving the properties for each molecular
"""
from pikachu.general import draw_smiles
from rdkit.Chem import AllChem
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
class reaction ():
    def __init__(self,substrates = "", products = ""):
        self.substrates = substrates
        self.products = products
    def get_reaction_site(self):
        bit_fingerprint = np.zeros((0,),
                                   dtype=int)  # (one dimention, 0 is number of rows)


        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius,
                                                                  num_bits)
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint)



def main():
    unittest.main()

if __name__ == "__main__":
    main()

