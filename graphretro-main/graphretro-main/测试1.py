from rdkit import Chem
from rdkit.Chem import  Draw
# from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults

from typing import List, Dict, Tuple, Set
from collections import namedtuple, deque
ReactionInfo = namedtuple("ReactionInfo", ['rxn_smi', 'core', 'core_edits', 'lg_edits', 'attach_atoms', 'rxn_class'])

from seq_graph_retro.utils.chem import apply_edits_to_mol, get_mol, get_sub_mol
from seq_graph_retro.molgraph.mol_features import BOND_FLOAT_TO_TYPE
from seq_graph_retro.molgraph import MultiElement
from seq_graph_retro.utils.parse import get_reaction_info, extract_leaving_groups,get_reaction_core,get_bond_info




rxn_smi = '[NH2:3][c:4]1[cH:5][cH:6][c:7]([O:8][c:9]2[cH:10][cH:11][n:12][c:13]3[nH:14][cH:15][cH:16][c:17]23)[c:18]([F:19])[cH:20]1.[O:1]=[C:2]([C:21]([F:22])([F:23])[F:24])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31]>>[O:1]=[C:2]([NH:3][c:4]1[cH:5][cH:6][c:7]([O:8][c:9]2[cH:10][cH:11][n:12][c:13]3[nH:14][cH:15][cH:16][c:17]23)[c:18]([F:19])[cH:20]1)[C:21]([F:22])([F:23])[F:24]'
r, p = rxn_smi.split(">>")
reac_mol = get_mol(r)
prod_mol = get_mol(p)

mols = [
    reac_mol,
    prod_mol
]
img=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,500),legends=['' for x in mols])
img.save("image_r_p.png")


rxn_core, core_edits = get_reaction_core(r, p, kekulize=True, use_h_labels=True)  # 获取反应中心

prod_bonds = get_bond_info(prod_mol)
p_amap_idx = {atom.GetAtomMapNum(): atom.GetIdx() for atom in prod_mol.GetAtoms()}

reaction_info = get_reaction_info(rxn_smi, kekulize=True,
                                  use_h_labels=True,
                                  rxn_class=int(5))




print(rxn_core,core_edits)
print('prod_bonds',prod_bonds)
print('p_amap_idx',p_amap_idx)
print('rxn',reaction_info[0])
print('反应中心',reaction_info[1])
print('核心编辑',reaction_info[2])
print('离去编辑',reaction_info[3])
print('附加原子',reaction_info[4])
print('反应类型',reaction_info[5])

############################################### 拆分产物p
fragments = apply_edits_to_mol(prod_mol, core_edits)
smiles = Chem.MolToSmiles(fragments)
r1, r2 = smiles.split(".")
smls = [r1,
        r2]
mols = []
for a in  smls:
    a = Chem.MolFromSmiles(a)
    mols.append(a)
img=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,500),legends=['' for x in mols])
img.save("chai_fen.png")
print(smiles)


