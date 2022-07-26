import pandas as pd
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, display, Image
from rdkit.Chem import Draw,DataStructs
from rdkit import Chem
import numpy as np
smile = "[1c:1]1(=[20O:20])[2nH:2][18c:18]([19NH2:19])[3n:3][4c:4]2[5c:5]1[6n+:6]([42CH3:42])[7cH:7][8n:8]2[9C@@H:9]1[10O:10][11C@H:11]([17CH2:17][16O:16][29P:29]([30O:30][31P:31]([32O:32][33P:33]([34O:34][41CH2:41][21C@H:21]2[22O:22][25C@@H:25]([26*:26])[24C@H:24]([27O:27][43CH3:43])[23C@@H:23]2[28O:28][53P:53]([54O:54][57CH2:57][44C@H:44]2[45O:45][48C@@H:48]([49*:49])[47C@H:47]([50OH:50])[46C@@H:46]2[51O:51][52*:52])(=[55O:55])[56O-:56])(=[37O:37])[40O-:40])(=[36O:36])[39O-:39])(=[35O:35])[38O-:38])[12C@@H:12]([13OH:13])[14C@H:14]1[15OH:15]"
mol = Chem.MolFromSmiles(smile)
Draw.ShowMol(mol, size=(600, 600))
radius = 3
atom_index = 50
atom_environment = rdmolops.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index)
atom_map = {}

submol = Chem.PathToSubmol(mol, atom_environment, atomMap=atom_map)
Chem.SanitizeMol(submol,
                 sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
rdmolops.SanitizeFlags.SANITIZE_NONE

Draw.ShowMol(submol, size=(600, 600))
bit_fingerprint_mol = np.zeros((0,),
                                   dtype=int)  # (one dimention, 0 is number of rows)

#returns an RDKit vector object.
substrate_molecular = submol
num_bits = 2048
morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(substrate_molecular, radius,
                                                          num_bits)

# We convert the RDKit vetor object to a numpy array.
DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint_mol)
print( sum(bit_fingerprint_mol))