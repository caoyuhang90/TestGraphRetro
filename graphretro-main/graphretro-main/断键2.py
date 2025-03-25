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
rxn_smi = ('[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31].[NH2:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1>>[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[NH:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1')
# rxn_smi = '[C](=[O])([C]([F])([F])[F])[O][C](=[O])[C]([F])([F])[F].[NH2][c]1[cH][cH][c]([O][c]2[cH][cH][n][c]3[nH][cH][cH][c]23)[c]([F])[cH]1>>[C](=[O])([C]([F])([F])[F])[NH][c]1[cH][cH][c]([O][c]2[cH][cH][n][c]3[nH][cH][cH][c]23)[c]([F])[cH]1'

r, p = rxn_smi.split(">>")
print("1111",p)
reac_mol = get_mol(r)
prod_mol = get_mol(p)
# print(len(prod_mol))
mols = [
    reac_mol,
    prod_mol
]
img=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,500),legends=['' for x in mols])
img.save("aaa_duanjian/image_r_p_0.png")



rxn_smi_new = remove_amap_not_in_product(rxn_smi)  # 修正编号
rxn_smi_new, _ = remap_rxn_smi(rxn_smi_new)  # 重新映射，让原子一一对应

r, p = rxn_smi_new.split(">>")
print("r:",r)
print('p:',p)


reac_mol = get_mol(r)
prod_mol = get_mol(p)

mols = [
    reac_mol,
    prod_mol
]
img=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,500),legends=['' for x in mols])
img.save("aaa_duanjian/image_r_p_1.png")
#
rxn_core, core_edits = get_reaction_core(r, p, kekulize=True, use_h_labels=True)
print('rxn_core',rxn_core)
print('core_edits',core_edits)




############################## 拆分测试
num1 = 8
num2 = 9
num3 = 1
num4 = 0
core_edits = [f'{num1}:{num2}:{num3}:{num4}']    #  a1:a2:b1:b2  --->   （产物原子编号）-起始：终止：原来的键类型：后来的键类型
fragments = apply_edits_to_mol(prod_mol, core_edits)  # 拆分--产物的分子表示+中心编辑
fragments = Chem.MolToSmiles(fragments)
print(fragments)
fragment_list = fragments.split('.')
molecules = [Chem.MolFromSmiles(fragment) for fragment in fragment_list]
img = Draw.MolsToGridImage(molecules, molsPerRow=1, subImgSize=(500, 500))
img.save("aaa_duanjian/chai_r_p_2.png")




