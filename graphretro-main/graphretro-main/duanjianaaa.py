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
    pcanon = canonicalize(p)  # 清除氢原子和其原子映射
    pmol = Chem.MolFromSmiles(pcanon)
    [atom.SetAtomMapNum(atom.GetIdx()+1) for atom in pmol.GetAtoms()]  # 设置原子映射编号
    p = Chem.MolToSmiles(pmol)
    return p

# rxn_smi = '[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31].[NH2:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1>>[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[NH:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1'
# rxn_smi = ('[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[O:27][C:26](=[O:25])[C:28]([F:29])([F:30])[F:31].[NH2:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1>>[C:1](=[O:2])([C:3]([F:4])([F:5])[F:6])[NH:7][c:8]1[cH:9][cH:10][c:11]([O:12][c:13]2[cH:14][cH:15][n:16][c:17]3[nH:18][cH:19][cH:20][c:21]23)[c:22]([F:23])[cH:24]1')
# rxn_smi = '[C](=[O])([C]([F])([F])[F])[O][C](=[O])[C]([F])([F])[F].[NH2][c]1[cH][cH][c]([O][c]2[cH][cH][n][c]3[nH][cH][cH][c]23)[c]([F])[cH]1>>[C](=[O])([C]([F])([F])[F])[NH][c]1[cH][cH][c]([O][c]2[cH][cH][n][c]3[nH][cH][cH][c]23)[c]([F])[cH]1'



######################################################  包装函数  用于 _init__.py 的 238 行附近
def new_edits_split(n1,n2=0,n3=-1,n4=0):
    num1 = n1  # 5  6
    num2 = n2
    num3 = 1
    num4 = 0

    num5 = n3  # 8 9
    num6 = n4
    num7 = 1
    num8 = 0


    core_edits = [f'{num1}:{num2}:{num3}:{num4}']    #  a1:a2:b1:b2  --->   （产物原子编号）-起始：终止：原来的键类型：后来的键类型
    core_edits_1 = [f'{num5}:{num6}:{num7}:{num8}']

    fragments = apply_edits_to_mol(prod_mol, core_edits)  # 拆分--产物的分子表示+中心编辑
    fragments = Chem.MolToSmiles(fragments)
    fragment_list = []
    fragment_list = fragments.split('.')



    # 双键操作
    for fra in fragment_list:
        fra = Chem.MolFromSmiles(fra)
        franum_list = []
        print('ok')
        for atom in fra.GetAtoms():
            atomnum = atom.GetAtomMapNum()
            franum_list.append(atomnum)
        if num5 in franum_list :  # 判断是在哪个子图里面
            # fragment_list = fragment_list.remove(Chem.MolToSmiles(fra))

            # 删除 符合的第一次拆分的子图
            str1_mol = Chem.MolToSmiles(fra)
            str1 = str(str1_mol)
            fragment_list.remove(str1)
            fragments_1 = apply_edits_to_mol(fra, core_edits_1)  # 拆分--产物的分子表示+中心编辑
            fragments_1 = Chem.MolToSmiles(fragments_1).split('.')
            fragment_list = fragment_list + fragments_1


    combined_smiles = '.'.join(fragment_list)
    print('combined_smiles',combined_smiles)
    return Chem.MolFromSmiles(combined_smiles)











############################# 拆分测试   拆不了【环】
def edits_split(prod_mol,n1,n2=0,n3=-1,n4=0):
    prod_mol = canonicalize_prod(prod_mol)
    prod_mol = get_mol(prod_mol)

    num1 = n1  # 5  6
    num2 = n2
    num3 = 1
    num4 = 0

    num5 = n3  # 8 9
    num6 = n4
    num7 = 1
    num8 = 0


    core_edits = [f'{num1}:{num2}:{num3}:{num4}']    #  a1:a2:b1:b2  --->   （产物原子编号）-起始：终止：原来的键类型：后来的键类型
    core_edits_1 = [f'{num5}:{num6}:{num7}:{num8}']

    fragments = apply_edits_to_mol(prod_mol, core_edits)  # 拆分--产物的分子表示+中心编辑
    fragments = Chem.MolToSmiles(fragments)
    fragment_list = []
    fragment_list = fragments.split('.')
    print('fragment_list',fragment_list)

    fragment_list_new = []
    # 双键操作
    for fra in fragment_list:
        fra = Chem.MolFromSmiles(fra)
        franum_list = []
        print('ok')
        for atom in fra.GetAtoms():
            atomnum = atom.GetAtomMapNum()
            franum_list.append(atomnum)
        if num5 in franum_list :  # 判断是在哪个子图里面
            # fragment_list = fragment_list.remove(Chem.MolToSmiles(fra))

            # 删除 符合的第一次拆分的子图
            str1_mol = Chem.MolToSmiles(fra)
            str1 = str(str1_mol)

            str1 = canonicalize(str1)
            print('str1', str1)

            # fragment_list_new = [ canonicalize(x)  for x in fragment_list]
            # fragment_list_new.remove(str1)

            fragments_1 = apply_edits_to_mol(fra, core_edits_1)  # 拆分--产物的分子表示+中心编辑
            print('fragments_1:',Chem.MolToSmiles(fragments_1))
            fragments_1 = Chem.MolToSmiles(fragments_1).split('.')
            print(fragments_1)
            fragment_list_new.extend(fragments_1)
            # print('fragment_list_aaa:',fragment_list)
        else:
            fra = Chem.MolToSmiles(fra)
            fragment_list_new.append(fra)

    highlight_atoms_list = []
    molecules = [Chem.MolFromSmiles(fragment) for fragment in fragment_list_new]
    for mol in molecules:
        highlight_atoms = []
        for atom in mol.GetAtoms():
            print(atom.GetAtomMapNum())
            if atom.GetAtomMapNum() in [n1, n2,n3,n4]:  # 目标映射编号
                highlight_atoms.append(atom.GetIdx())
        highlight_atoms_list.append(highlight_atoms)


    img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(800, 500),highlightAtomLists=highlight_atoms_list)
    img.save("aaa_duanjian/chai_r_p_3_2.png")


# [H][C@]12[C@@]3(C[C@@H](C4[C@](C=CC(C=4)=O)([C@]3([C@H](C[C@@]1([C@]([C@H](C)C2)(C(=O)SCF)OC(CC)=O)C)O)F)C)F)[H]
# rxn_smi = ('CCNC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1/C=C/[C@H](CCC2=CC=CC=C2)O)O)O')
# rxn_smi = ('O=C(C)C1=CC=C(CN2C=C(C=N2)NC(C2=C(C3C=CC=CC=3)OC(COC)=N2)=O)O1')
# rxn_smi = '[H][C@]12[C@@]3(C[C@@H](C4[C@](C=CC(C=4)=O)([C@]3([C@H](C[C@@]1([C@]([C@H](C)C2)(C(=O)SCF)OC(CC)=O)C)O)F)C)F)[H]'  # 5 6 15 16
rxn_smi = 'CCNC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1/C=C/[C@H](CCC2=CC=CC=C2)O)O)O'






# edits_split(rxn_smi,9,10,14,15)
# new = new_edits_split(5,6,15,16)

print('okduanjian')



