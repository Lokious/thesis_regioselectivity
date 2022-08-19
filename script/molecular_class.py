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
            mol = Chem.MolFromSmiles(r"{}".format(input))
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
        This function is to create MorganFingerprint from mol
        
        :param substrate_molecular: mol, substrate of the methylation reaction
        :param num_bits:int, the length of fingerprint
        :param radius:int use for calculate fingerprint (“radius” from an atom is measured by the number of bonds that separates two atoms)
        :return:bit_fingerprint_mol: np.array, MorganFingerprint
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

    def create_fingerprint_atom(self,  substrate_molecular: Chem.Mol,
        atom_object: Chem.Atom, num_bits: int = 2048,radius: int = 3) -> np.array:
        """
        Create MorganFingerprint for substructure in the molecule of input atom

        :param substrate_molecular: mol, substrate of the methylation reaction
        :param atom_object: Atom, one atom from substrate_molecular
        :param num_bits:int, the length of fingerprint
        :param radius:int use for calculate fingerprint
        :return: bit_fingerprint_atom: np.array, MorganFingerprint
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
        raise RuntimeError("function `perform_reaction()` is deprecated")
        
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
        :param product_smiles:sting,the SMILES of products
            substrate_smiles:sting,the SMILES of substrates
        :return:
        """
        #remove AtomMapNum and Isotope, avoiding affect the similarity
        product_mols=[]
        substrate_mols=[]
        for item in product_smiles.split("."):
            product_mol = Chem.MolFromSmiles(r"{}".format(item))
            # sanitize molecular
            Chem.SanitizeMol(
                product_mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            rdmolops.SanitizeFlags.SANITIZE_NONE
            for atom in product_mol.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)
            product_mols.append(product_mol)
        for item in substrate_smiles.split("."):
            substrate_mol = Chem.MolFromSmiles(r"{}".format(item))
            # sanitize molecular
            Chem.SanitizeMol(
                substrate_mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            rdmolops.SanitizeFlags.SANITIZE_NONE
            for atom in substrate_mol.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)
            substrate_mols.append(substrate_mol)
        #Draw.ShowMol(product_mol,(600,600))

        #get editable mol for removing methylation group
        mol = Chem.EditableMol(product_mol)
        atom_index_list = []
        similiarity = 0
        #key is mol object value is possiable atom index
        possiable_dictionary={}
        for i,product_mol in enumerate(product_mols):
            for atom in product_mol.GetAtoms():
                #do not consider C-methylation here
                # if atom.GetSymbol() !="C":
                #     neighbours = atom.GetNeighbors()
                #     for neighbour_atom in neighbours:
                #find -CH3 group
                if ((atom.GetSymbol()=="C") and (atom.GetDegree()==1)):

                    #try to remove this -CH3 group
                    mol = Chem.EditableMol(product_mol)
                    mol.RemoveAtom(atom.GetIdx())
                    #save possiable atoms to a list, if remove one do not get 100% similarity substrate
                    #the list will be use later

                    # Draw.ShowMol(mol.GetMol(), (600, 600))
                    # Draw.ShowMol(substrate_mol, (600, 600))
                    #compare similarity between substrate and molecule after remove -CH3
                    for j,substrate_mol in enumerate(substrate_mols):
                        if (len(product_mol.GetAtoms())-len(substrate_mol.GetAtoms()))>4:
                            continue
                        else:

                            try:
                                possiable_dictionary[(i, j)].append(
                                    atom.GetIdx())
                            except:
                                possiable_dictionary[(i, j)]=[]
                                possiable_dictionary[(i,j)].append(atom.GetIdx())

                            similiarity = DataStructs.FingerprintSimilarity(
                                Chem.RDKFingerprint(mol.GetMol()),
                                Chem.RDKFingerprint(substrate_mol))
                            print(similiarity)
                            #if similarity equals 1, it matches the substrate,
                            if similiarity==1:
                                #rest isotope
                                neighbour_atom=(product_mol.GetAtomWithIdx(atom.GetIdx())).GetNeighbors()[0]
                                atom_index_list.append(
                                    neighbour_atom.GetSymbol() + ":"+str(
                                        neighbour_atom.GetIdx()))
                                for atom1 in product_mol.GetAtoms():
                                    atom1.SetIsotope(atom1.GetIdx())
                                #reset the mol and remove the methlation group, keep the isotope
                                mol=Chem.EditableMol(product_mol)
                                mol.RemoveAtom(atom.GetIdx())
                                mol_remove_methylation = mol.GetMol()
                                #Draw.ShowMol(mol_remove_methylation, (600, 600))
                                print("check")
                                return product_mol,Chem.MolToSmiles(mol_remove_methylation),atom_index_list,"Pass_Check"
                            else:
                                #mol=mol_save_status
                                continue
        #remove multiple methylation group
        print(possiable_dictionary)
        from itertools import combinations
        for key in possiable_dictionary.keys():
            #possiable list is the
            possiable_list=possiable_dictionary[key]
            index_combination_list=list(combinations(possiable_list, 2))
            index_combination_list += list(combinations(possiable_list, 3))
            print(index_combination_list)
            print(key)
            product_mol=product_mols[key[0]]
            for items in index_combination_list:
                mol=Chem.EditableMol(product_mol)
                remove_atom=sorted(list(items),reverse=True)
                #remove from larger index
                for i in remove_atom:
                    mol.RemoveAtom(i)
                #Draw.ShowMol(mol.GetMol(), (600, 600))
                #Draw.ShowMol(product_mol,(600,600))
                #Draw.ShowMol(substrate_mols[key[1]],(600,600))
                similiarity = DataStructs.FingerprintSimilarity(
                    Chem.RDKFingerprint(mol.GetMol()),
                    #get the substrate_mol through index which is saved in the key of dictionary
                    Chem.RDKFingerprint(substrate_mols[key[1]]))
                print(similiarity)
                if similiarity ==1:
                    for i in remove_atom:
                        atom = (product_mol.GetAtomWithIdx(i)).GetNeighbors()[0]
                        atom_index_list.append(
                            atom.GetSymbol() + ":"+str(
                                atom.GetIdx()))
                    print(atom_index_list)
                        # rest isotope
                    for atom1 in product_mol.GetAtoms():
                        atom1.SetIsotope(atom1.GetIdx())
                    # reset the mol and remove the methlation group, keep the isotope
                    mol = Chem.EditableMol(product_mol)
                    mol.RemoveAtom(atom.GetIdx())
                    mol_remove_methylation = mol.GetMol()
                    #Draw.ShowMol(mol.GetMol(), (600, 600))
                    print("check")
                    return product_mol, Chem.MolToSmiles(
                        mol_remove_methylation), atom_index_list, "Pass_Check"

        mol_remove_methylation=mol.GetMol()

        #print(similiarity)
        #raise RuntimeError("function `get_reaction_sites()` is deprecated")

        return product_mol, Chem.MolToSmiles(
                        mol_remove_methylation), atom_index_list, "unCheck"
        
    def get_reactant_atom(self):
        """
        This function only consider the large substrate and product, it assume
        the smaller molecular as Methyl donor, H2O or coenzyme etc.

        :return: list of atom objects which is the regioselectivity site
        """
        raise RuntimeError("function `get_reactant_atom()` is deprecated")
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
    raise RuntimeError("function `main_substrate()` is deprecated")
    sim_dictionary = {}
    mol_object = Molecule()
    reaction_object = Reaction()

    for i, mol1_smile in enumerate(subs):

        try:
            mol1 = Chem.MolFromSmiles(r"{}".format(mol1_smile))
            for atom in mol1.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)

        except:
            return None

        for j,mol2_smile in enumerate(pros):

            # mol2 = Chem.MolFromSmiles(r"{}".format(mol2_smile))
            # Draw.ShowMol(mol2)

            try:

                mol2 = Chem.MolFromSmiles(r"{}".format(mol2_smile))
                for atom in mol2.GetAtoms():
                    atom.SetAtomMapNum(0)
                    atom.SetIsotope(0)

            except:
                return None
            from rdkit.Chem.Draw import SimilarityMaps
            # target_mol_simi_fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
            #     mol2, mol1, SimilarityMaps.GetMorganFingerprint)
            #target_mol_simi_fig.show()
            sim_dictionary[(i, j)] = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2))


    similarity_list_top = list(
        sorted(
            sim_dictionary.items(),
            reverse=True,
            key=lambda item: item[1]
        )
    )

    for key in similarity_list_top:
        i = key[0][0]
        j= key[0][1]

        # print(len(pros[j]))
        # print(len(subs[i]))
        # the documentation below is not quite clear to me

        #this function assumed the similarity of fingerprint betwween molecular before methylation
        #and after methylation should be higher than this molecular with other molecular
        #the methyl donor molecular will become smaller after reaction
        mol_pro=Chem.MolFromSmiles(pros[j])
        atom_numb_pro=len(mol_pro.GetAtoms())
        print(atom_numb_pro)
        #Draw.ShowMol(mol_pro)
        mol_sub=Chem.MolFromSmiles(subs[j])
        atom_numb_sub=len(mol_sub.GetAtoms())
        #Draw.ShowMol(mol_sub)
        print(atom_numb_sub)
        if (len(pros[j]) > len(subs[i])) and (pros[j] != 0) and (subs[i] != 0 ):
            mol1=Chem.MolFromSmiles(r"{}".format(subs[i]))
            mol2 = Chem.MolFromSmiles(r"{}".format(pros[j]))
            # Draw.ShowMol(mol2)
            # Draw.ShowMol(mol1)
            return [subs[i], pros[j]]
        else:
            continue
                
                
class Testreaction_class(unittest.TestCase):
    def test0_get_reaction_sites(self):
        """
        Test if get_reaction_sites works for single reactant
        """
        reaction = Reaction()
        substrates="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1OH:7])[4cH:4][2cH:2]1"
        products="c1[1c:1]([6CH2:8][7CH:9]=[8C:10]([9CH3:11])[10CH3:12])[3c:3]([OH:6])[5c:5]([1O:7][11CH3:13])[4cH:4][2cH:2]1"

        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(list_methylsite[0], "O7:1")

    def test1_get_reactant_atom(self):
        """
        Test for get_reaction_sites works for multiple methylation site
        """
        reaction = Reaction()
        substrates = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n:5][3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30OH:85])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([51CH3:116])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30O:85][50CH3:115])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"
        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(list_methylsite,['O:40', 'N:8'])
        self.assertEqual(check,"pass check")

    def test2_get_reactant_atom(self):
        """
        Exception situation, when valency change in reaction

        """
        reaction = Reaction()
        substrates = r"O[As:1]([CH3:2])[1OH:3].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"[As](=[O:1])([1O-:2])([CH3:3])[1CH3:4].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"
        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertNotEqual(list_methylsite,['As:1'])
        self.assertEqual(check,'unCheck')

def main():
    #unittest.main()
    reaction= Reaction()
    substrates = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n:5][3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30OH:85])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
    products = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([51CH3:116])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30O:85][50CH3:115])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"
    reaction.get_reaction_sites(products,substrates)

    # mian_list = main_substrate(subs=substrates, pros=products)
    # print(mian_list)
    # reaction = Reaction()
    # df = pd.read_csv("../data/seq_smiles_all_test_codefor_finding_methylsite.csv",header=0,index_col=0)
    # df.dropna(inplace=True)
    # df_save=copy.deepcopy(df)
    # df_save["remove_methylation"]=pd.DataFrame(
    #         len(df_save.index) * [0]).astype('object')
    # df_save["list_methylsite"]=pd.DataFrame(
    #         len(df_save.index) * [0]).astype('object')
    # df_save["check"]=pd.DataFrame(
    #         len(df_save.index) * [0]).astype('object')
    # for index in df.index:
    #     print(index)
    #     main_sub_smile = df.loc[index,"sub_smiles"]
    #     main_pro_smile = df.loc[index,"pro_smiles"]
    #     pro_mol,remove_methyl_smile,list_methylsite,check=reaction.get_reaction_sites(main_pro_smile,main_sub_smile)
    #     if check=="pass check":
    #         Draw.MolToFile(Chem.MolFromSmiles(remove_methyl_smile), "../pro_fig/{}_remove.png".format(index))
    #     # Draw.MolToFile(pro_mol,
    #     #                "../pro_fig/{}_pro.png".format(index))
    #     df_save.loc[index,"remove_methylation"]=remove_methyl_smile
    #     df_save.loc[index,"list_methylsite"]=",".join(list_methylsite)
    #     df_save.loc[index,"check"]=check
    # print(df_save)
    # df_save.to_csv("test_finding_methylation_site.csv")

if __name__ == "__main__":
    main()
