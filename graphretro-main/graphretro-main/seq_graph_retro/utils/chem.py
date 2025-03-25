from rdkit import Chem
from typing import Iterable, List

# Redfined from seq_graph_retro/molgraph/mol_features.py
BOND_TYPES = [None, Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, \
    Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]
BOND_FLOAT_TO_TYPE = {
    0.0: BOND_TYPES[0],
    1.0: BOND_TYPES[1],
    2.0: BOND_TYPES[2],
    3.0: BOND_TYPES[3],
    1.5: BOND_TYPES[4],
}

def get_mol(smiles: str, kekulize: bool = False) -> Chem.Mol:
    """SMILES string to Mol.

    Parameters
    ----------
    smiles: str,
        SMILES string for molecule
    kekulize: bool,
        Whether to kekulize the molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and kekulize: # kekulize：是否清除芳香环信息
        Chem.Kekulize(mol)
    return mol

def apply_edits_to_mol(mol: Chem.Mol, edits: Iterable[str]) -> Chem.Mol:
    """Apply edits to molecular graph.

    Parameters
    ----------
    mol: Chem.Mol,
        RDKit mol object
    edits: Iterable[str],
        Iterable of edits to apply. An edit is structured as a1:a2:b1:b2, where
        a1, a2 are atom maps of participating atoms and b1, b2 are previous and
        new bond orders. When  a2 = 0, we update the hydrogen count.

        a1：第一个参与原子的映射编号（原子的索引）。
        a2：第二个参与原子的映射编号。如果 a2 为 0，则表示要更新第一个原子的氢原子数量。
        b1：键的原始顺序（类型），用于检查或验证。
        b2：新的键的顺序（类型），表示将要设置的新值。

    """
    new_mol = Chem.RWMol(mol)  # RWMol 允许进行原子和键的添加、删除和修改
    amap = {atom.GetAtomMapNum(): atom.GetIdx() for atom in new_mol.GetAtoms()}  # amap { 原子编号：原子ID }

    # Keep track of aromatic nitrogens, might cause explicit hydrogen issues
    aromatic_nitrogen_idx = set()  # 存储芳香氮原子的索引。
    aromatic_carbonyl_adj_to_aromatic_nH = {}  # 存储与芳香氮原子相邻的碳酰基原子。
    aromatic_carbondeg3_adj_to_aromatic_nH0 = {}  # 存储与芳香氮原子相邻的三价碳原子。
    for a in new_mol.GetAtoms():  # 新分子中的每个原子
        if a.GetIsAromatic() and a.GetSymbol() == 'N':  # 判断原子是否是芳香性原子 且  是否为N原子
            aromatic_nitrogen_idx.add(a.GetIdx())   #  记录芳香氮原子索引
            for nbr in a.GetNeighbors():   # 遍历氮原子的邻近原子，检查是否是碳原子、是否芳香、是否有双键、是否有三个键。
                nbr_is_carbon = (nbr.GetSymbol() == 'C')
                nbr_is_aromatic = nbr.GetIsAromatic()
                nbr_has_double_bond = any(b.GetBondTypeAsDouble() == 2 for b in nbr.GetBonds())
                nbr_has_3_bonds = (len(nbr.GetBonds()) == 3)

                # 根据氮原子的氢原子数量和邻近碳原子的性质，将相关索引存入字典。
                if (a.GetNumExplicitHs() ==1 and nbr_is_carbon and nbr_is_aromatic
                    and nbr_has_double_bond):
                    aromatic_carbonyl_adj_to_aromatic_nH[nbr.GetIdx()] = a.GetIdx()
                elif (a.GetNumExplicitHs() == 0 and nbr_is_carbon and nbr_is_aromatic
                    and nbr_has_3_bonds):
                    aromatic_carbondeg3_adj_to_aromatic_nH0[nbr.GetIdx()] = a.GetIdx()
        else:
            a.SetNumExplicitHs(0)  # 非芳香氮原子，设置其显性氢原子数量为 0
    new_mol.UpdatePropertyCache()  # 更新分子属性缓存  ？？

    # Apply the edits as predicted
    for edit in edits:
        x, y, prev_bo, new_bo = edit.split(":")  # 起始：结束：原来的键类型：后来的键类型
        x, y = int(x), int(y)
        new_bo = float(new_bo)

        if y == 0:
            continue

        print('x, y, prev_bo, new_bo',x, y, prev_bo, new_bo)
        # x和y是编号，索引其对于的ID号，针对ID号操作
        bond = new_mol.GetBondBetweenAtoms(amap[x],amap[y])  # 获取x与y之间的键类型
        a1 = new_mol.GetAtomWithIdx(amap[x])
        a2 = new_mol.GetAtomWithIdx(amap[y])

        if bond is not None:
            new_mol.RemoveBond(amap[x],amap[y])  # 移除原来的键类型

            # Are we losing a bond on an aromatic nitrogen?
            ''' 判断键类型
            
            单键（SINGLE）：返回 1.0
            双键（DOUBLE）：返回 2.0
            三键（TRIPLE）：返回 3.0
            芳香键（AROMATIC）：返回 1.5
            其他键类型：根据具体情况返回相应的数值
            
            '''

            if bond.GetBondTypeAsDouble() == 1.0:  # 单键
                if amap[x] in aromatic_nitrogen_idx:
                    if a1.GetTotalNumHs() == 0:
                        a1.SetNumExplicitHs(1)
                    elif a1.GetFormalCharge() == 1:
                        a1.SetFormalCharge(0)
                elif amap[y] in aromatic_nitrogen_idx:
                    if a2.GetTotalNumHs() == 0:
                        a2.SetNumExplicitHs(1)
                    elif a2.GetFormalCharge() == 1:
                        a2.SetFormalCharge(0)

            # Are we losing a c=O bond on an aromatic ring? If so, remove H from adjacent nH if appropriate
            if bond.GetBondTypeAsDouble() == 2.0:  # 双键
                if amap[x] in aromatic_carbonyl_adj_to_aromatic_nH:
                    new_mol.GetAtomWithIdx(aromatic_carbonyl_adj_to_aromatic_nH[amap[x]]).SetNumExplicitHs(0)
                elif amap[y] in aromatic_carbonyl_adj_to_aromatic_nH:
                    new_mol.GetAtomWithIdx(aromatic_carbonyl_adj_to_aromatic_nH[amap[y]]).SetNumExplicitHs(0)

        if new_bo > 0:  # 添加新键
            new_mol.AddBond(amap[x],amap[y],BOND_FLOAT_TO_TYPE[new_bo])  # BOND_FLOAT_TO_TYPE[new_bo]：将new_bo浮点数转化为BondType类型

            # Special alkylation case?
            if new_bo == 1:
                if amap[x] in aromatic_nitrogen_idx:
                    if a1.GetTotalNumHs() == 1:
                        a1.SetNumExplicitHs(0)
                    else:
                        a1.SetFormalCharge(1)
                elif amap[y] in aromatic_nitrogen_idx:
                    if a2.GetTotalNumHs() == 1:
                        a2.SetNumExplicitHs(0)
                    else:
                        a2.SetFormalCharge(1)

            # Are we getting a c=O bond on an aromatic ring? If so, add H to adjacent nH0 if appropriate
            if new_bo == 2:
                if amap[x] in aromatic_carbondeg3_adj_to_aromatic_nH0:
                    new_mol.GetAtomWithIdx(aromatic_carbondeg3_adj_to_aromatic_nH0[amap[x]]).SetNumExplicitHs(1)
                elif amap[y] in aromatic_carbondeg3_adj_to_aromatic_nH0:
                    new_mol.GetAtomWithIdx(aromatic_carbondeg3_adj_to_aromatic_nH0[amap[y]]).SetNumExplicitHs(1)

    pred_mol = new_mol.GetMol()  # 可编辑的分子对象 new_mol 获取最终的分子对象 pred_mol

    # Clear formal charges to make molecules valid
    # Note: because S and P (among others) can change valence, be more flexible
    for atom in pred_mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1: # exclude negatively-charged azide
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            if bond_vals <= 3:
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'N' and atom.GetFormalCharge() == -1: # handle negatively-charged azide addition
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            if bond_vals == 3 and any([nbr.GetSymbol() == 'N' for nbr in atom.GetNeighbors()]):
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'N':
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            if bond_vals == 4 and not atom.GetIsAromatic(): # and atom.IsInRingSize(5)):
                atom.SetFormalCharge(1)
        elif atom.GetSymbol() == 'C' and atom.GetFormalCharge() != 0:
            atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'O' and atom.GetFormalCharge() != 0:
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]) + atom.GetNumExplicitHs()
            if bond_vals == 2:
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() in ['Cl', 'Br', 'I', 'F'] and atom.GetFormalCharge() != 0:
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            if bond_vals == 1:
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'S' and atom.GetFormalCharge() != 0:
            bond_vals = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            if bond_vals in [2, 4, 6]:
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'P': # quartenary phosphorous should be pos. charge with 0 H
            bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            if sum(bond_vals) == 4 and len(bond_vals) == 4:
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(0)
            elif sum(bond_vals) == 3 and len(bond_vals) == 3: # make sure neutral
                atom.SetFormalCharge(0)
        elif atom.GetSymbol() == 'B': # quartenary boron should be neg. charge with 0 H
            bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            if sum(bond_vals) == 4 and len(bond_vals) == 4:
                atom.SetFormalCharge(-1)
                atom.SetNumExplicitHs(0)
        elif atom.GetSymbol() in ['Mg', 'Zn']:
            bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            if sum(bond_vals) == 1 and len(bond_vals) == 1:
                atom.SetFormalCharge(1)
        elif atom.GetSymbol() == 'Si':
            bond_vals = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            if sum(bond_vals) == len(bond_vals):
                atom.SetNumExplicitHs(max(0, 4 - len(bond_vals)))

    return pred_mol

def get_sub_mol(mol: Chem.Mol, sub_atoms: List[int]) -> Chem.Mol:
    """Extract subgraph from molecular graph.

    Parameters
    ----------
    mol: Chem.Mol,
        RDKit mol object,
    sub_atoms: List[int],
        List of atom indices in the subgraph.
    """
    new_mol = Chem.RWMol()
    atom_map = {}
    for idx in sub_atoms:
        atom = mol.GetAtomWithIdx(idx)
        atom_map[idx] = new_mol.AddAtom(atom)

    sub_atoms = set(sub_atoms)
    for idx in sub_atoms:
        a = mol.GetAtomWithIdx(idx)
        for b in a.GetNeighbors():
            if b.GetIdx() not in sub_atoms: continue
            bond = mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            bt = bond.GetBondType()
            if a.GetIdx() < b.GetIdx(): #each bond is enumerated twice
                new_mol.AddBond(atom_map[a.GetIdx()], atom_map[b.GetIdx()], bt)

    return new_mol.GetMol()

def get_sub_mol_stereo(mol: Chem.Mol, sub_atoms: List[int]) -> Chem.Mol:
    """Extract subgraph from molecular graph, while preserving stereochemistry.

    Parameters
    ----------
    mol: Chem.Mol,
        RDKit mol object,
    sub_atoms: List[int],
        List of atom indices in the subgraph.
    """
    # This version retains stereochemistry, as opposed to the other version
    new_mol = Chem.RWMol(Chem.Mol(mol))
    atoms_to_remove = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in sub_atoms:
            atoms_to_remove.append(atom.GetIdx())

    for atom_idx in sorted(atoms_to_remove, reverse=True):
        new_mol.RemoveAtom(atom_idx)

    return new_mol.GetMol()
