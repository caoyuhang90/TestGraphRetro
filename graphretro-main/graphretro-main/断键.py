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
from seq_graph_retro.utils.edit_mol import canonicalize

from data_process.canonicalize_prod import remove_amap_not_in_product,remap_rxn_smi


def canonicalize_prod(p):
    pcanon = canonicalize(p)  # 清楚氢原子和其原子映射
    pmol = Chem.MolFromSmiles(pcanon)
    [atom.SetAtomMapNum(atom.GetIdx()+1) for atom in pmol.GetAtoms()]  # 设置原子映射编号
    p = Chem.MolToSmiles(pmol)
    return p

# rxn_smi = '[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31].[NH2:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1>>[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[NH:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1'
rxn_smi = ('[NH2:3][c:4]1[cH:5][cH:6][c:7]([O:8][c:9]2[cH:10][cH:11][n:12][c:13]3[nH:14][cH:15][cH:16][c:17]23)[c:18]([F:19])[cH:20]1.[O:1]=[C:2]([C:21]([F:22])([F:23])[F:24])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31]>>'
           '[O:1]=[C:2]([NH:3][c:4]1[cH:5][cH:6][c:7]([O:8][c:9]2[cH:10][cH:11][n:12][c:13]3[nH:14][cH:15][cH:16][c:17]23)[c:18]([F:19])[cH:20]1)[C:21]([F:22])([F:23])[F:24]')

r, p = rxn_smi.split(">>")
p = canonicalize_prod(p)
r_can = canonicalize(r)

print(p,r_can)


print(p)
print('r',r)
print('r_can',r_can)

reac_mol = get_mol(r_can)
prod_mol = get_mol(p)

mols = [
    reac_mol,
    prod_mol
]
img=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,500),legends=['' for x in mols])
img.save("aaa_duanjian/image_r_p.png")

rxn_core, core_edits = get_reaction_core(r, p, kekulize=True, use_h_labels=True)
print(rxn_core,core_edits)


