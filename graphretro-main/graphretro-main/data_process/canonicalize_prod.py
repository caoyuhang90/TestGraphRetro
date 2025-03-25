"""
Canonicalize the product SMILES, and then use substructure matching to infer
the correspondence to the original atom-mapped order. This correspondence is then
used to renumber the reactant atoms.
"""

from rdkit import Chem
import os
import argparse
import pandas as pd

#DATA_DIR = f"{os.environ['SEQ_GRAPH_RETRO']}/datasets/uspto-50k/"
DATA_DIR = r"D:\11111AAA_Retro_pro\GraphRetro(1)\graphretro-main\graphretro-main\datasets\uspto-50k"

# 把skiles标准化并编号
def canonicalize_prod(p):
    import copy
    p = copy.deepcopy(p)  # # 深拷贝输入，以避免修改原始数据
    p = canonicalize(p)  # 调用 canonicalize 函数对 SMILES 进行标准化
    p_mol = Chem.MolFromSmiles(p)   # # 将标准化后的 SMILES 转换为分子对象
    for atom in p_mol.GetAtoms():  # 遍历分子中的每个原子
        atom.SetAtomMapNum(atom.GetIdx() + 1)  # 编号
    p = Chem.MolToSmiles(p_mol)   # 将分子对象转换回 SMILES 表示
    return p

def remove_amap_not_in_product(rxn_smi):
    """
    Corrects the atom map numbers of atoms only in reactants. 
    This correction helps avoid the issue of duplicate atom mapping
    after the canonicalization step.

    目的是修正反应物中的原子映射编号，以避免在标准化步骤后出现重复的原子映射。
    """
    r, p = rxn_smi.split(">>")  # 将反应物和产物分开

    pmol = Chem.MolFromSmiles(p)   # 将产物的 SMILES 转换为分子对象
    pmol_amaps = set([atom.GetAtomMapNum() for atom in pmol.GetAtoms()])  # # 获取产物中所有原子的映射编号
    max_amap = max(pmol_amaps) #Atoms only in reactants are labelled starting with max_amap  找最大的编号

    rmol  = Chem.MolFromSmiles(r)

    for atom in rmol.GetAtoms():  # 遍历反应物中的每个原子
        amap_num = atom.GetAtomMapNum()   # 获取当前原子的映射编号
        if amap_num not in pmol_amaps:  # 如果该映射编号不在产物的原子映射编号中
            atom.SetAtomMapNum(max_amap+1)  #  # 设置新的映射编号
            max_amap += 1

    r_updated = Chem.MolToSmiles(rmol)   # 将更新后的反应物分子对象转换为 SMILES
    rxn_smi_updated = r_updated + ">>" + p  # 重新组合反应物和产物的 SMILES
    return rxn_smi_updated

# 标准化：移除氢原子，移除所有的原子编号
def canonicalize(smiles):
    try:
        tmp = Chem.MolFromSmiles(smiles)  # 尝试将 SMILES 转换为分子对象
    except:
        print('no mol', flush=True)
        return smiles
    if tmp is None:
        return smiles
    tmp = Chem.RemoveHs(tmp)  # # 移除分子中的氢原子
    [a.ClearProp('molAtomMapNumber') for a in tmp.GetAtoms()]  # 清除原子映射编号属性
    return Chem.MolToSmiles(tmp)

def infer_correspondence(p):
    orig_mol = Chem.MolFromSmiles(p)  # 原始分子
    canon_mol = Chem.MolFromSmiles(canonicalize_prod(p))  # 标准化分子
    matches = list(canon_mol.GetSubstructMatches(orig_mol))[0]  # 原始分子和标准化分子进行匹配   #这里分子太长子结构匹配会进入死循环

    # 创建原始分子的字典{原子的索引：原子的编号}
    idx_amap = {atom.GetIdx(): atom.GetAtomMapNum() for atom in orig_mol.GetAtoms()}

    correspondence = {}
    for idx, match_idx in enumerate(matches):
        match_anum = canon_mol.GetAtomWithIdx(match_idx).GetAtomMapNum()
        old_anum = idx_amap[idx]
        correspondence[old_anum] = match_anum
    return correspondence


def remap_rxn_smi(rxn_smi):
    r, p = rxn_smi.split(">>")
    canon_mol = Chem.MolFromSmiles(canonicalize_prod(p))
    correspondence = infer_correspondence(p)

    rmol = Chem.MolFromSmiles(r)

    # 把前体的编号和产物的编号统一起来，让产物原子追踪到前体原子编号，让他们一一对应
    # 确保反应物中每个原子的映射编号与标准化产物中的原子映射编号相对应。保持原子之间的对应关系。
    for atom in rmol.GetAtoms():
        atomnum = atom.GetAtomMapNum()
        if atomnum in correspondence:
            newatomnum = correspondence[atomnum]
            atom.SetAtomMapNum(newatomnum)

    rmol = Chem.MolFromSmiles(Chem.MolToSmiles(rmol))
    rxn_smi_new = Chem.MolToSmiles(rmol) + ">>" + Chem.MolToSmiles(canon_mol)
    return rxn_smi_new, correspondence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default=DATA_DIR, help="Directory where data is located.")
    parser.add_argument("--filename", required=True, help="File with reactions to canonicalize")
    args = parser.parse_args()

    new_file = f"canonicalized_{args.filename}"
    df = pd.read_csv(f"{args.data_dir}/{args.filename}")
    print(f"Processing file of size: {len(df)}")

    new_dict = {'id': [], 'class': [], 'reactants>reagents>production': []}
    for idx in range(len(df)):
        element = df.loc[idx]  # 读取 行 数据
        uspto_id, class_id, rxn_smi = element['id'], element['class'], element['reactants>reagents>production']
        
        rxn_smi_new = remove_amap_not_in_product(rxn_smi)  # 修正编号
        rxn_smi_new, _ = remap_rxn_smi(rxn_smi_new)  # 重新映射，让原子一一对应
        new_dict['id'].append(uspto_id)
        new_dict['class'].append(class_id)
        new_dict['reactants>reagents>production'].append(rxn_smi_new)

    new_df = pd.DataFrame.from_dict(new_dict)
    new_df.to_csv(f"{args.data_dir}/{new_file}", index=False)

if __name__ == "__main__":
    main()
