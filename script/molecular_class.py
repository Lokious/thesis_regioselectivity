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
from rdkit.Chem import MACCSkeys
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

    MACCSkey_dictionary = {
        1: ('?', 0),  # ISOTOPE
        # 2:('[#104,#105,#106,#107,#106,#109,#110,#111,#112]',0),  # atomic num >103 Not complete
        2: ('[#104]', 0),
        # limit the above def'n since the RDKit only accepts up to #104
        3: ('[#32,#33,#34,#50,#51,#52,#82,#83,#84]', 0),
        # Group IVa,Va,VIa Rows 4-6
        4: ('[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]', 0),  # actinide
        5: ('[Sc,Ti,Y,Zr,Hf]', 0),  # Group IIIB,IVB (Sc...)
        6: ('[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]', 0),  # Lanthanide
        7: ('[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]', 0),  # Group VB,VIB,VIIB
        8: ('[!#6;!#1]1~*~*~*~1', 0),  # QAAA@1
        9: ('[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]', 0),  # Group VIII (Fe...)
        10: ('[Be,Mg,Ca,Sr,Ba,Ra]', 0),  # Group IIa (Alkaline earth)
        11: ('*1~*~*~*~1', 0),  # 4M Ring
        12: ('[Cu,Zn,Ag,Cd,Au,Hg]', 0),  # Group IB,IIB (Cu..)
        13: ('[#8]~[#7](~[#6])~[#6]', 0),  # ON(C)C
        14: ('[#16]-[#16]', 0),  # S-S
        15: ('[#8]~[#6](~[#8])~[#8]', 0),  # OC(O)O
        16: ('[!#6;!#1]1~*~*~1', 0),  # QAA@1
        17: ('[#6]#[#6]', 0),  # CTC
        18: ('[#5,#13,#31,#49,#81]', 0),  # Group IIIA (B...)
        19: ('*1~*~*~*~*~*~*~1', 0),  # 7M Ring
        20: ('[#14]', 0),  # Si
        21: ('[#6]=[#6](~[!#6;!#1])~[!#6;!#1]', 0),  # C=C(Q)Q
        22: ('*1~*~*~1', 0),  # 3M Ring
        23: ('[#7]~[#6](~[#8])~[#8]', 0),  # NC(O)O
        24: ('[#7]-[#8]', 0),  # N-O
        25: ('[#7]~[#6](~[#7])~[#7]', 0),  # NC(N)N
        26: ('[#6]=;@[#6](@*)@*', 0),  # C$=C($A)$A
        27: ('[I]', 0),  # I
        28: ('[!#6;!#1]~[CH2]~[!#6;!#1]', 0),  # QCH2Q
        29: ('[#15]', 0),  # P
        30: ('[#6]~[!#6;!#1](~[#6])(~[#6])~*', 0),  # CQ(C)(C)A
        31: ('[!#6;!#1]~[F,Cl,Br,I]', 0),  # QX
        32: ('[#6]~[#16]~[#7]', 0),  # CSN
        33: ('[#7]~[#16]', 0),  # NS
        34: ('[CH2]=*', 0),  # CH2=A
        35: ('[Li,Na,K,Rb,Cs,Fr]', 0),  # Group IA (Alkali Metal)
        36: ('[#16R]', 0),  # S Heterocycle
        37: ('[#7]~[#6](~[#8])~[#7]', 0),  # NC(O)N
        38: ('[#7]~[#6](~[#6])~[#7]', 0),  # NC(C)N
        39: ('[#8]~[#16](~[#8])~[#8]', 0),  # OS(O)O
        40: ('[#16]-[#8]', 0),  # S-O
        41: ('[#6]#[#7]', 0),  # CTN
        42: ('F', 0),  # F
        43: ('[!#6;!#1;!H0]~*~[!#6;!#1;!H0]', 0),  # QHAQH
        44: ('[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]', 0),
        # OTHER
        45: ('[#6]=[#6]~[#7]', 0),  # C=CN
        46: ('Br', 0),  # BR
        47: ('[#16]~*~[#7]', 0),  # SAN
        48: ('[#8]~[!#6;!#1](~[#8])(~[#8])', 0),  # OQ(O)O
        49: ('[!+0]', 0),  # CHARGE
        50: ('[#6]=[#6](~[#6])~[#6]', 0),  # C=C(C)C
        51: ('[#6]~[#16]~[#8]', 0),  # CSO
        52: ('[#7]~[#7]', 0),  # NN
        53: ('[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]', 0),  # QHAAAQH
        54: ('[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]', 0),  # QHAAQH
        55: ('[#8]~[#16]~[#8]', 0),  # OSO
        56: ('[#8]~[#7](~[#8])~[#6]', 0),  # ON(O)C
        57: ('[#8R]', 0),  # O Heterocycle
        58: ('[!#6;!#1]~[#16]~[!#6;!#1]', 0),  # QSQ
        59: ('[#16]!:*:*', 0),  # Snot%A%A
        60: ('[#16]=[#8]', 0),  # S=O
        61: ('*~[#16](~*)~*', 0),  # AS(A)A
        62: ('*@*!@*@*', 0),  # A$!A$A
        63: ('[#7]=[#8]', 0),  # N=O
        64: ('*@*!@[#16]', 0),  # A$A!S
        65: ('c:n', 0),  # C%N
        66: ('[#6]~[#6](~[#6])(~[#6])~*', 0),  # CC(C)(C)A
        67: ('[!#6;!#1]~[#16]', 0),  # QS
        68: ('[!#6;!#1;!H0]~[!#6;!#1;!H0]', 0),  # QHQH (&...) SPEC Incomplete
        69: ('[!#6;!#1]~[!#6;!#1;!H0]', 0),  # QQH
        70: ('[!#6;!#1]~[#7]~[!#6;!#1]', 0),  # QNQ
        71: ('[#7]~[#8]', 0),  # NO
        72: ('[#8]~*~*~[#8]', 0),  # OAAO
        73: ('[#16]=*', 0),  # S=A
        74: ('[CH3]~*~[CH3]', 0),  # CH3ACH3
        75: ('*!@[#7]@*', 0),  # A!N$A
        76: ('[#6]=[#6](~*)~*', 0),  # C=C(A)A
        77: ('[#7]~*~[#7]', 0),  # NAN
        78: ('[#6]=[#7]', 0),  # C=N
        79: ('[#7]~*~*~[#7]', 0),  # NAAN
        80: ('[#7]~*~*~*~[#7]', 0),  # NAAAN
        81: ('[#16]~*(~*)~*', 0),  # SA(A)A
        82: ('*~[CH2]~[!#6;!#1;!H0]', 0),  # ACH2QH
        83: ('[!#6;!#1]1~*~*~*~*~1', 0),  # QAAAA@1
        84: ('[NH2]', 0),  # NH2
        85: ('[#6]~[#7](~[#6])~[#6]', 0),  # CN(C)C
        86: ('[C;H2,H3][!#6;!#1][C;H2,H3]', 0),  # CH2QCH2
        87: ('[F,Cl,Br,I]!@*@*', 0),  # X!A$A
        88: ('[#16]', 0),  # S
        89: ('[#8]~*~*~*~[#8]', 0),  # OAAAO
        90:
            (
            '[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]',
            0),  # QHAACH2A
        91:
            (
            '[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]',
            0),  # QHAAACH2A
        92: ('[#8]~[#6](~[#7])~[#6]', 0),  # OC(N)C
        93: ('[!#6;!#1]~[CH3]', 0),  # QCH3
        94: ('[!#6;!#1]~[#7]', 0),  # QN
        95: ('[#7]~*~*~[#8]', 0),  # NAAO
        96: ('*1~*~*~*~*~1', 0),  # 5 M ring
        97: ('[#7]~*~*~*~[#8]', 0),  # NAAAO
        98: ('[!#6;!#1]1~*~*~*~*~*~1', 0),  # QAAAAA@1
        99: ('[#6]=[#6]', 0),  # C=C
        100: ('*~[CH2]~[#7]', 0),  # ACH2N
        101:
            (
            '[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]',
            0),  # 8M Ring or larger. This only handles up to ring sizes of 14
        102: ('[!#6;!#1]~[#8]', 0),  # QO
        103: ('Cl', 0),  # CL
        104: ('[!#6;!#1;!H0]~*~[CH2]~*', 0),  # QHACH2A
        105: ('*@*(@*)@*', 0),  # A$A($A)$A
        106: ('[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]', 0),  # QA(Q)Q
        107: ('[F,Cl,Br,I]~*(~*)~*', 0),  # XA(A)A
        108: ('[CH3]~*~*~*~[CH2]~*', 0),  # CH3AAACH2A
        109: ('*~[CH2]~[#8]', 0),  # ACH2O
        110: ('[#7]~[#6]~[#8]', 0),  # NCO
        111: ('[#7]~*~[CH2]~*', 0),  # NACH2A
        112: ('*~*(~*)(~*)~*', 0),  # AA(A)(A)A
        113: ('[#8]!:*:*', 0),  # Onot%A%A
        114: ('[CH3]~[CH2]~*', 0),  # CH3CH2A
        115: ('[CH3]~*~[CH2]~*', 0),  # CH3ACH2A
        116: ('[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]', 0),  # CH3AACH2A
        117: ('[#7]~*~[#8]', 0),  # NAO
        118: ('[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]', 1),  # ACH2CH2A > 1
        119: ('[#7]=*', 0),  # N=A
        120: ('[!#6;R]', 1),  # Heterocyclic atom > 1 (&...) Spec Incomplete
        121: ('[#7;R]', 0),  # N Heterocycle
        122: ('*~[#7](~*)~*', 0),  # AN(A)A
        123: ('[#8]~[#6]~[#8]', 0),  # OCO
        124: ('[!#6;!#1]~[!#6;!#1]', 0),  # QQ
        125: ('?', 0),  # Aromatic Ring > 1
        126: ('*!@[#8]!@*', 0),  # A!O!A
        127: ('*@*!@[#8]', 1),  # A$A!O > 1 (&...) Spec Incomplete
        128:
            (
            '[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]',
            0),  # ACH2AAACH2A
        129: (
        '[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]',
        0),  # ACH2AACH2A
        130: ('[!#6;!#1]~[!#6;!#1]', 1),  # QQ > 1 (&...)  Spec Incomplete
        131: ('[!#6;!#1;!H0]', 1),  # QH > 1
        132: ('[#8]~*~[CH2]~*', 0),  # OACH2A
        133: ('*@*!@[#7]', 0),  # A$A!N
        134: ('[F,Cl,Br,I]', 0),  # X (HALOGEN)
        135: ('[#7]!:*:*', 0),  # Nnot%A%A
        136: ('[#8]=*', 1),  # O=A>1
        137: ('[!C;!c;R]', 0),  # Heterocycle
        138: ('[!#6;!#1]~[CH2]~*', 1),  # QCH2A>1 (&...) Spec Incomplete
        139: ('[O;!H0]', 0),  # OH
        140: ('[#8]', 3),  # O > 3 (&...) Spec Incomplete
        141: ('[CH3]', 2),  # CH3 > 2  (&...) Spec Incomplete
        142: ('[#7]', 1),  # N > 1
        143: ('*@*!@[#8]', 0),  # A$A!O
        144: ('*!:*:*!:*', 0),  # Anot%A%Anot%A
        145: ('*1~*~*~*~*~*~1', 1),  # 6M ring > 1
        146: ('[#8]', 2),  # O > 2
        147: ('[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]', 0),  # ACH2CH2A
        148: ('*~[!#6;!#1](~*)~*', 0),  # AQ(A)A
        149: ('[C;H3,H4]', 1),  # CH3 > 1
        150: ('*!@*@*!@*', 0),  # A!A$A!A
        151: ('[#7;!H0]', 0),  # NH
        152: ('[#8]~[#6](~[#6])~[#6]', 0),  # OC(C)C
        153: ('[!#6;!#1]~[CH2]~*', 0),  # QCH2A
        154: ('[#6]=[#8]', 0),  # C=O
        155: ('*!@[CH2]!@*', 0),  # A!CH2!A
        156: ('[#7]~*(~*)~*', 0),  # NA(A)A
        157: ('[#6]-[#8]', 0),  # C-O
        158: ('[#6]-[#7]', 0),  # C-N
        159: ('[#8]', 1),  # O>1
        160: ('[C;H3,H4]', 0),  # CH3
        161: ('[#7]', 0),  # N
        162: ('a', 0),  # Aromatic
        163: ('*1~*~*~*~*~*~1', 0),  # 6M Ring
        164: ('[#8]', 0),  # O
        165: ('[R]', 0),  # Ring
        166: ('?', 0),  # Fragments  FIX: this can't be done in SMARTS
    }
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
        bi = {}
        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(
            substrate_molecular, 
            radius,
            num_bits,bitInfo=bi
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
        bi = {}
        #returns an RDKit vector object.
        morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(
            submol, 
            radius,
         num_bits,bitInfo=bi
        )
        
        # We convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_atom)

        return bit_fingerprint_atom

    def create_MACCSkey_mol(self, substrate_molecular: Chem.Mol) -> np.array:

        """
        This function is to get MACCSkeys for substrate molecule

        https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/MACCSkeys.py
        :param substrate_molecular:  mol, substrate of the methylation reaction
        :return:
        """
        # sanitize molecular
        Chem.SanitizeMol(
            substrate_molecular,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
        rdmolops.SanitizeFlags.SANITIZE_NONE
        # initialize a numpy array for molecular fingerprint

        bit_fingerprint_mol = np.zeros(
            (0,),
            dtype=int  # (one dimention, 0 is number of rows)
        )

        MACCSkey_bit_vector = MACCSkeys.GenMACCSKeys(substrate_molecular)

        DataStructs.ConvertToNumpyArray(MACCSkey_bit_vector,
                                        bit_fingerprint_mol)
        print(bit_fingerprint_mol)

        return bit_fingerprint_mol

    def create_MACCSkey_atom(self, substrate_molecular: Chem.Mol,
                             atom_object: Chem.Atom,
                             radius: int = 3) -> np.array:
        atom_index = atom_object.GetIdx()
        atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(
            substrate_molecular,
            radius,
            atom_index
        )
        atom_map = {}
        submol = Chem.PathToSubmol(substrate_molecular, atom_environment,
                                   atomMap=atom_map)
        # sanitize molecule
        Chem.SanitizeMol(submol,
                         sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        rdmolops.SanitizeFlags.SANITIZE_NONE
        # inisilize a numpy array for molecular fingerprint
        bit_fingerprint_atom = np.zeros(
            (0,),
            dtype=int  # (one dimention, 0 is number of rows)
        )
        # returns an RDKit vector object(MACCSKey).
        MACCSkey_bit_vector = MACCSkeys.GenMACCSKeys(submol)

        # We convert the RDKit vetor object to a numpy array.
        DataStructs.ConvertToNumpyArray(MACCSkey_bit_vector,
                                        bit_fingerprint_atom)
        print(bit_fingerprint_atom)
        return bit_fingerprint_atom


class Reaction:
    def __init__(self, substrates="", products="", rxn_object=None):
        self.substrates = substrates
        self.products = products
        self.rxn_object = rxn_object
        self.mol_product = None
        self.mol_substrate = None


    def get_reaction_sites(self,product_smiles:str="",substrate_smiles:str=""):
        """
        finding methylation atoms by removing -CH3 group and compare the similarity between substrate and product

        :param product_smiles:sting,the SMILES of products, split by '.'
            substrate_smiles:sting,the SMILES of substrates, split by '.'
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
                #find -CH3 group
                if ((atom.GetSymbol()=="C") and (atom.GetDegree()==1)):

                    #try to remove this -CH3 group
                    mol = Chem.EditableMol(product_mol)
                    mol.RemoveAtom(atom.GetIdx())
                    #compare similarity between substrate and molecule after remove -CH3
                    for j,substrate_mol in enumerate(substrate_mols):
                        if (len(product_mol.GetAtoms())-len(substrate_mol.GetAtoms()))>4:
                            continue
                        else:
                            # save possiable atoms to a list, if remove one do not get 100% similarity substrate
                            # the list will be use later
                            try:
                                possiable_dictionary[(i, j)].append(
                                    atom.GetIdx())
                            except:
                                possiable_dictionary[(i, j)]=[]
                                possiable_dictionary[(i,j)].append(atom.GetIdx())

                            similiarity = DataStructs.FingerprintSimilarity(
                                Chem.RDKFingerprint(mol.GetMol()),
                                Chem.RDKFingerprint(substrate_mol))
                            #print(similiarity)
                            #if similarity equals 1, it matches the substrate,
                            if similiarity==1:
                                #rest isotope
                                neighbour_atom=(product_mol.GetAtomWithIdx(atom.GetIdx())).GetNeighbors()[0]
                                atom_index_list.append(
                                    neighbour_atom.GetSymbol() + ":"+str(
                                        neighbour_atom.GetIdx()))
                                atom_index=neighbour_atom.GetIdx()
                                # check radical electrons
                                print("symbol:{},index{}, Explicit valence: {},total H:{},TotalDegree{}".format(neighbour_atom.GetSymbol(),neighbour_atom.GetIdx(),
                                      neighbour_atom.GetExplicitValence(),
                                      neighbour_atom.GetTotalNumHs(),neighbour_atom.GetTotalDegree()))
                                for atom1 in product_mol.GetAtoms():
                                    atom1.SetIsotope(atom1.GetIdx())
                                #reset the mol and remove the methlation group, keep the isotope
                                mol=Chem.EditableMol(product_mol)
                                mol.RemoveAtom(atom.GetIdx())

                                mol_remove_methylation = mol.GetMol()
                                # Add H after remove methyl group
                                for atom in mol_remove_methylation.GetAtoms():
                                    if atom.GetIsotope() == atom_index:
                                        # Add H after remove methyl group
                                        num_H=atom.GetNumExplicitHs()
                                        #print("num_H:{}".format(atom.GetNumExplicitHs()))
                                        try:
                                            atom.SetNumExplicitHs((atom.GetExplicitValence() - atom.GetTotalDegree()))
                                            similiarity = DataStructs.FingerprintSimilarity(
                                                Chem.RDKFingerprint(
                                                    mol_remove_methylation),
                                                # get the substrate_mol through index which is saved in the key of dictionary
                                                Chem.RDKFingerprint(substrate_mol))
                                            print(similiarity)
                                            assert similiarity == 1
                                            Chem.Kekulize(mol_remove_methylation)
                                        except:
                                            #print("num_H:{}".format(atom.GetNumExplicitHs()))
                                            atom.SetNumExplicitHs(num_H)
                                            #print("num_H:{}".format(atom.GetNumExplicitHs()))
                                            print("Kekulize error after add H, won't add H to the atom")
                                Draw.ShowMol(product_mol)
                                Draw.ShowMol(mol_remove_methylation, (800, 800))
                                #Draw.ShowMol(substrate_mol)
                                print("check")
                                return product_mol,Chem.MolToSmiles(mol_remove_methylation),atom_index_list,"Pass_Check"
                            else:
                                #mol=mol_save_status
                                continue
        #remove multiple methylation group
        print(possiable_dictionary)
        from itertools import combinations
        for key in possiable_dictionary.keys():
            #possiable list is the atom index of the carbon in the methyl group
            possiable_list=possiable_dictionary[key]
            index_combination_list=list(combinations(possiable_list, 2))
            index_combination_list += list(combinations(possiable_list, 3))
            # print(index_combination_list)
            # print(key)
            product_mol=product_mols[key[0]]
            for items in index_combination_list:
                mol=Chem.EditableMol(product_mol)
                remove_atom=sorted(list(items),reverse=True)
                #remove from larger index
                for i in remove_atom:
                    mol.RemoveAtom(i)
                similiarity = DataStructs.FingerprintSimilarity(
                    Chem.RDKFingerprint(mol.GetMol()),
                    #get the substrate_mol through index which is saved in the key of dictionary
                    Chem.RDKFingerprint(substrate_mols[key[1]]))

                print(similiarity)
                if similiarity ==1:
                    # Draw.ShowMol(product_mol, (800, 800))
                    # Draw.ShowMol(mol.GetMol(), (800, 800))
                    # Draw.ShowMol(substrate_mols[key[1]], (800, 800))

                    # reset the mol and remove the methlation group, keep the isotope
                    for atom1 in product_mol.GetAtoms():
                        # rest isotope
                        atom1.SetIsotope(atom1.GetIdx())
                    #convert mol to editable mol
                    mol = Chem.EditableMol(product_mol)
                    for i in remove_atom:
                        atom = (product_mol.GetAtomWithIdx(i)).GetNeighbors()[0]
                        atom_index_list.append(
                            atom.GetSymbol() + ":"+str(
                                atom.GetIdx()))
                        mol.RemoveAtom(i)
                        # Add H after remove methyl group
                        try:
                            print("Total degree{}".format(atom.GetTotalDegree()))
                            atom.SetNumExplicitHs((atom.GetExplicitValence() - atom.GetTotalDegree()))
                            Chem.Kekulize(mol.GetMol())
                        except:
                            print(
                            "Kekulize error after add H, won't add H to the atom")
                    #print(atom_index_list)
                    mol_remove_methylation = mol.GetMol()
                    # Draw.ShowMol(mol_remove_methylation, (800, 800))
                    print("check")
                    if DataStructs.FingerprintSimilarity(
                    Chem.RDKFingerprint(mol_remove_methylation),
                    #get the substrate_mol through index which is saved in the key of dictionary
                    Chem.RDKFingerprint(substrate_mols[key[1]]))==1:
                        return product_mol, Chem.MolToSmiles(
                            mol_remove_methylation), atom_index_list, "Pass_Check"
                    else:
                        raise ValueError
        mol_remove_methylation=mol.GetMol()

        return product_mol, Chem.MolToSmiles(
                        mol_remove_methylation), atom_index_list, "unCheck"

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
            mol1 = Chem.MolFromSmiles(r"{}".format(subs[i]))
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
        self.assertEqual(list_methylsite[0], "O:10")
        self.assertEqual(check, "Pass_Check")

    def test1_get_reaction_sites(self):
        """
        Test for get_reaction_sites works for multiple methylation site
        """
        reaction = Reaction()
        substrates = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n:5][3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30OH:85])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([51CH3:116])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:20]([5O:21][1P:22]([6O:23][2P:24](=[9O:27])([12O-:30])[15O:39][14CH2:38][12C@H:35]2[13O:33][10C@@H:32]([5n:31]3[32c:81]4[11n:80][31cH:79][10n:78][30c:77]([9NH2:76])[33c:82]4[12n:83][34cH:84]3)[11C@H:34]([30O:85][50CH3:115])[13C@@H:36]2[14O:37][3P:40]([16O:41][19CH2:51][17C@H:48]2[19O:46][15C@@H:44]([6n:45]3[42c:102]4[16n:101][41cH:100][15n:99][40c:98]([14NH2:97])[43c:103]4[17n:104][44cH:105]3)[16C@H:47]([31OH:86])[18C@@H:49]2[20O:50][4P:52]([21O:53][22CH2:60][21C@H:59]2[24O:58][20C@@H:57]([7n:56]3[35c:89](=[34O:92])[13nH:91][38c:94](=[35O:96])[37c:93]([39CH3:95])[36cH:90]3)[24C@H:63]([32OH:87])[23C@@H:61]2[25O:62][5P:64]([26O:65][29CH2:75][27C@H:72]2[29O:70][25C@@H:69]([8n:68]3[47c:111]4[20n:110][46cH:109][19n:108][45c:107]([18NH2:106])[48c:112]4[21n:113][49cH:114]3)[26C@H:71]([33OH:88])[28C@@H:73]2[*:74])(=[27O:66])[28O-:67])(=[22O:54])[23O-:55])(=[17O:42])[18O-:43])(=[8O:26])[11O-:29])(=[7O:25])[10O-:28])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"

        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(list_methylsite[0],'O:40')
        self.assertEqual(list_methylsite[1], 'N:8')
        self.assertEqual(check,"Pass_Check")


    def test2_get_reaction_sites(self):
        """
        Exception situation, when valency change in reaction

        """
        reaction = Reaction()
        substrates = r"O[As:1]([CH3:2])[1OH:3].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"[As](=[O:1])([1O-:2])([CH3:3])[1CH3:4].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"
        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(list_methylsite,[])
        self.assertEqual(check,'unCheck')
    def test3_get_reaction_sites(self):
        """
        This is test for add H after remove the methyl group
        :return:
        """
        reaction = Reaction()
        substrates = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([15CH3:42])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:29]([8O:30][1P:31]([9O:32][2P:33]([10O:34][14CH2:41][10C@H:20]2[5O:21][13C@@H:24]([*:25])[12C@H:23]([6OH:26])[11C@@H:22]2[7O:27][1*:28])(=[13O:37])[16O-:40])(=[12O:36])[15O-:39])(=[11O:35])[14O-:38])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([15CH3:42])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:29]([8O:30][1P:31]([9O:32][2P:33]([10O:34][14CH2:41][10C@H:20]2[5O:21][13C@@H:24]([*:25])[12C@H:23]([6O:26][16CH3:43])[11C@@H:22]2[7O:27][1*:28])(=[13O:37])[16O-:40])(=[12O:36])[15O-:39])(=[11O:35])[14O-:38])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"

        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(check,"Pass_Check")
        mol = Chem.MolFromSmiles(remove_methyl_smile)
        for atom in mol.GetAtoms():
            if atom.GetIsotope()==29:
                num_electrons=atom.GetNumRadicalElectrons()
                explicitValence=atom.GetExplicitValence()
                self.assertEqual(num_electrons,0)
                self.assertEqual(explicitValence,2)
    def test4_get_reaction_sites(self):
        """
        Test for get_reaction_sites works for not adding H
        """
        reaction = Reaction()
        substrates = r"O[C@H]1C[C@@H](O[C@@H]1COP([O-])([O-])=O)n1ccc(=O)[nH]c1=O"
        products = r"Cc1cn([C@H]2C[C@H](O)[C@@H](COP([O-])([O-])=O)O2)c(=O)[nH]c1=O"
        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        self.assertEqual(list_methylsite[0],'C:1')
        self.assertEqual(check,"Pass_Check")
        print(remove_methyl_smile)
    def test5_get_reaction_sites(self):
        print("TEST5")
        reaction = Reaction()
        substrates = r"c1[1cH:1][3c:5]([2NH2:18])[1n:4][2c:3](=[1O:10])[n:2]1[4C@@H:6]1[O:7][5C@H:8]([7CH2:15][4O:14][P:13](=[2O:11])([3O-:12])[1*:21])[6C@@H:9]([6O:19][*:20])[8C@H:16]1[5OH:17].C[S+:1]([1CH2:2][2CH2:3][3C@H:4]([NH3+:5])[4C:6]([O-:7])=[1O:8])[5CH2:9][6C@H:10]1[2O:11][7C@@H:12]([1n:17]2[10cH:18][2n:19][11c:20]3[12c:21]2[3n:22][13cH:23][4n:24][14c:25]3[5NH2:26])[8C@H:13]([3OH:14])[9C@@H:15]1[4OH:16]"
        products = r"c1[1cH:1][3c:5]([2NH:16][9CH3:22])[1n:4][2c:3](=[1O:10])[n:2]1[4C@@H:6]1[O:7][5C@H:8]([7CH2:15][4O:14][P:13](=[2O:11])([3O-:12])[*:20])[6C@@H:9]([6O:19][1*:21])[8C@H:17]1[5OH:18].[C@H]1([2OH:6])[1C@@H:1]([1OH:5])[2C@H:2]([n:7]2[4c:8]3[6c:10]([2n:12][5cH:9]2)[7c:13]([4NH2:16])[3n:15][8cH:14][1n:11]3)[O:3][3C@@H:4]1[9CH2:17][S:18][10CH2:19][11CH2:20][12C@@H:21]([13C:24]([3O-:22])=[4O:23])[5NH3+:25]"
        pro_mol,remove_methyl_smile,list_methylsite,check = reaction.get_reaction_sites(products,substrates)
        print(remove_methyl_smile)
        mol1=Chem.MolFromSmiles(remove_methyl_smile)
        Draw.ShowMol(mol1,(800,800))
        print(list_methylsite)
    def test_6_create_MACCSkey_mol(self):
        molecule=Molecule()
        substrates = r"O[C@H]1C[C@@H](O[C@@H]1COP([O-])([O-])=O)n1ccc(=O)[nH]c1=O"
        mol = Chem.MolFromSmiles(substrates)
        mol_fp=molecule.create_MACCSkey_mol(mol)
        print(len(mol_fp))
        self.assertEqual(mol_fp, 167)
    def test_7_create_MACCSkey_atom(self):
        molecule=Molecule()
        substrates = r"O[C@H]1C[C@@H](O[C@@H]1COP([O-])([O-])=O)n1ccc(=O)[nH]c1=O"
        mol = Chem.MolFromSmiles(substrates)
        for atom in mol.GetAtoms():
            atom_fp=molecule.create_MACCSkey_atom(mol,atom,3)
            #print(len(mol_fp))
            self.assertEqual(atom_fp,167)
def main():
    unittest.main()
    """
    reaction = Reaction()
    substrates = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([15CH3:42])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:29]([8O:30][1P:31]([9O:32][2P:33]([10O:34][14CH2:41][10C@H:20]2[5O:21][13C@@H:24]([*:25])[12C@H:23]([6OH:26])[11C@@H:22]2[7O:27][1*:28])(=[13O:37])[16O-:40])(=[12O:36])[15O-:39])(=[11O:35])[14O-:38])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14]"
    products = r"c1(=[4O:19])[nH:1][9c:17]([4NH2:18])[1n:2][1c:3]2[2c:4]1[2n+:5]([15CH3:42])[3cH:6][3n:7]2[4C@@H:8]1[O:9][5C@H:10]([8CH2:16][3O:15][P:29]([8O:30][1P:31]([9O:32][2P:33]([10O:34][14CH2:41][10C@H:20]2[5O:21][13C@@H:24]([*:25])[12C@H:23]([6O:26][16CH3:43])[11C@@H:22]2[7O:27][1*:28])(=[13O:37])[16O-:40])(=[12O:36])[15O-:39])(=[11O:35])[14O-:38])[6C@@H:11]([1OH:12])[7C@H:13]1[2OH:14]"
    #convert to mol object
    product_mol=Chem.MolFromSmiles(products)
    # clear the mapnumber and isotope
    for atom in product_mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)

    # set isotope
    for atom1 in product_mol.GetAtoms():
        atom1.SetIsotope(atom1.GetIdx())
    # draw mol
    Draw.ShowMol(product_mol,(600,600))
    # remove the methyl group
    mol = Chem.EditableMol(product_mol)

    mol.RemoveAtom(30)
    # convert to mol
    mol_remove=mol.GetMol()
    #Chem.SanitizeMol(mol_remove)

    print(mol_remove.GetAtomWithIdx(29).GetNumRadicalElectrons(),mol_remove.GetAtomWithIdx(29).GetTotalDegree(),mol_remove.GetAtomWithIdx(29).GetExplicitValence(),mol_remove.GetAtomWithIdx(29).GetImplicitValence(),mol_remove.GetAtomWithIdx(29).GetTotalValence())

    mol_remove.GetAtomWithIdx(29).SetNumExplicitHs((mol_remove.GetAtomWithIdx(29).GetExplicitValence()-mol_remove.GetAtomWithIdx(29).GetTotalDegree()))
    Draw.ShowMol(mol_remove,(600,600))
    print("symbol:{},index{}, Explicit valence: {},total H:{}".format(mol_remove.GetAtomWithIdx(29).GetSymbol(),mol_remove.GetAtomWithIdx(29).GetIsotope(),mol_remove.GetAtomWithIdx(29).GetExplicitValence(),mol_remove.GetAtomWithIdx(29).GetTotalNumHs()))
    #convert to smile
    substrate_smile = Chem.MolToSmiles(mol_remove)
    print(substrate_smile)
    #convert to mol object
    substrate_mol=Chem.MolFromSmiles(substrate_smile)
    #Chem.SanitizeMol(substrate_mol)
    Draw.ShowMol(substrate_mol,(600,600))
    for atom in substrate_mol.GetAtoms():
        if atom.GetIsotope()==29:

            print(atom.GetNumRadicalElectrons())
            print(
                "symbol:{},isotope(equal to before index){}, Explicit valence: {},total H:{}".format(
                    atom.GetSymbol(),
                    atom.GetIsotope(),
                    atom.GetExplicitValence(),
                    atom.GetTotalNumHs()))

    # pro_mol, remove_methyl_smile, list_methylsite, check = reaction.get_reaction_sites(
    #     products, substrates)
    # smile = remove_methyl_smile
    # print(smile)
    # print(substrate_smile)
    # print(smile==substrate_smile)
    # mol1= Chem.MolFromSmiles(smile)
    #
    # print(list_methylsite)
    # for atom in mol1.GetAtoms():
    #     if atom.GetIsotope()==29:
    #         index= atom.GetIdx()
    # print("symbol:{},isotope(equal to before index){}, Explicit valence: {},total H:{}".format(mol1.GetAtomWithIdx(index).GetSymbol(),mol1.GetAtomWithIdx(index).GetIsotope(),
    #       mol1.GetAtomWithIdx(index).GetExplicitValence(),
    #       mol1.GetAtomWithIdx(index).GetTotalNumHs()))
    # Draw.ShowMol(mol1,(600,600))
    """
if __name__ == "__main__":
    main()
