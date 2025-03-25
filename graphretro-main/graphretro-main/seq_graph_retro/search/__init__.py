import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem
import copy
import sys

from seq_graph_retro.utils.chem import apply_edits_to_mol
from seq_graph_retro.molgraph import RxnElement, MultiElement
from seq_graph_retro.data.collate_fns import pack_graph_feats, tensorize_bond_graphs
from seq_graph_retro.molgraph.mol_features import BOND_FLOATS
from duanjianaaa import new_edits_split

class BeamNode:

    def __init__(self, mol):
        self.mol = mol
        self._build_attrs()
        self._build_tensors()

    def _build_attrs(self):
        self.edit = None
        self.lg_groups = []
        self.lg_idx = 0
        self.prob = 0.0
        self.num_fragments = None
        self.node_complete = False

    def _build_tensors(self):
        self.prod_vecs = None
        self.frag_vecs = None
        self.prev_embed = None

    def add_edit(self, edit, edit_prob):
        self.edit = edit
        self.prob += edit_prob

    def add_lg(self, lg_group, lg_prob):
        if self.node_complete:
            print("All leaving groups added. Skipping for this node.")
            sys.stdout.flush()

        else:
            self.lg_groups.append(lg_group)
            self.lg_idx += 1
            self.prob += lg_prob
            if len(self.lg_groups) == self.num_fragments:
                self.node_complete = True

def copy_node(node):
    new_node = BeamNode(mol=Chem.Mol(node.mol))
    new_node.edit = node.edit
    new_node.lg_groups = copy.deepcopy(node.lg_groups)
    new_node.num_fragments = node.num_fragments
    new_node.prob = node.prob
    new_node.lg_idx = node.lg_idx
    new_node.node_complete = node.node_complete

    if node.prod_vecs is not None:
        new_node.prod_vecs = node.prod_vecs.clone()

    if node.frag_vecs is not None:
        new_node.frag_vecs = node.frag_vecs.clone()

    if node.prev_embed is not None:
        new_node.prev_embed = node.prev_embed.clone()

    return new_node

def check_nodes_complete(node_list):
    for node in node_list:
        if not node.node_complete:
            return False
    return True

class BeamSearch:

    def __init__(self, model, beam_width=1, max_edits=1):
        self.model = model
        self.beam_width = beam_width
        self.max_edits = max_edits

        if hasattr(self.model, 'encoder_name'):
            if self.model.encoder_name == 'GraphFeatEncoder':
                self.directed = True
            elif self.model.encoder_name == 'WLNEncoder':
                self.directed = False
            else:
                raise ValueError()
        else:
            if self.model.edit_net.encoder_name == 'GraphFeatEncoder':
                self.edit_directed = True
            elif self.model.edit_net.encoder_name == 'WLNEncoder':
                self.edit_directed = False
            else:
                raise ValueError()

            if self.model.lg_net.encoder_name == 'GraphFeatEncoder':
                self.lg_directed = True
            elif self.model.lg_net.encoder_name == 'WLNEncoder':
                self.lg_directed = False
            else:
                raise ValueError()

    def keep_topk_nodes(self, node_list):
        sorted_node_list = [copy_node(node) for node in node_list]
        sorted_node_list = sorted(sorted_node_list, key=lambda x: x.prob, reverse=True)
        if len(sorted_node_list) <= self.beam_width:
            return sorted_node_list

        sorted_node_list = sorted_node_list[:self.beam_width]
        return sorted_node_list

    def keep_topk_nodes_oneBroken(self, node_list,beam_number):
        sorted_node_list = [copy_node(node) for node in node_list]
        sorted_node_list = sorted(sorted_node_list, key=lambda x: x.prob, reverse=True)
        if len(sorted_node_list) <= beam_number:
            return sorted_node_list

        sorted_node_list = sorted_node_list[:beam_number]
        return sorted_node_list

    def keep_topk_nodes_new(self, node_list,second_num=5):
        indexed_node_list = [(index, copy_node(node)) for index, node in enumerate(node_list)]
        sorted_node_list = sorted(indexed_node_list, key=lambda x: x[1].prob, reverse=True)
        if len(sorted_node_list) <= second_num:
            return sorted_node_list

        sorted_node_list = sorted_node_list[:second_num]
        return sorted_node_list

    def get_topk_edits(self, prod_smi, rxn_class=None, **kwargs):
        mol = Chem.MolFromSmiles(prod_smi)
        node_list = [BeamNode(mol=Chem.Mol(mol)) for _ in range(self.beam_width)]

        use_rxn_class = False
        if rxn_class is not None:
            use_rxn_class = True

        prod_graph = RxnElement(mol=Chem.Mol(mol), rxn_class=rxn_class)
        directed = self.directed if hasattr(self, 'directed') else self.edit_directed

        if self.max_edits == 1:
            prod_tensors, prod_scopes = pack_graph_feats([prod_graph],
                                                         directed=directed,
                                                         return_graphs=False,
                                                         use_rxn_class=use_rxn_class)

            prod_tensors = self.model.to_device(prod_tensors)
            bg_inputs = tensorize_bond_graphs([prod_graph], directed=directed,
                                               use_rxn_class=use_rxn_class)
            bg_tensors, bg_scope = bg_inputs
            bg_tensors = self.model.to_device(bg_tensors)
            bg_inputs = (bg_tensors, bg_scope)
            prod_vecs, edit_logits, _ = self.model._compute_edit_logits(prod_tensors, prod_scopes,
                                                                        ha=None, bg_inputs=bg_inputs)
            edit_logits = edit_logits[0]

            edit_logits = F.log_softmax(edit_logits, dim=-1)
            k = min(len(edit_logits), self.beam_width)
            topk_vals, topk_idxs = torch.topk(edit_logits, k=k)

            for beam_idx, (topk_idx, val) in enumerate(zip(*(topk_idxs, topk_vals))):
                edit = self.get_edit_from_logits(mol=Chem.Mol(mol),
                                                 edit_logits=edit_logits,
                                                 idx=topk_idx, val=val)
                if not isinstance(edit, list):
                    edit = [edit]
                node_list[beam_idx].add_edit(edit=edit, edit_prob=val.item())
                if hasattr(self.model, 'encoder'):
                    node_list[beam_idx].prod_vecs = prod_vecs.clone()

            return [copy_node(node) for node in node_list]
        else:
            raise ValueError("Greater than 1 sequence length not supported yet.")

    def remove_invalid_nodes(self, prod_smi, node_list, rxn_class=None):
        new_list = []

        use_rxn_class = False
        if rxn_class is not None:
            use_rxn_class = True

        for node in node_list:
            try:
                mol = Chem.MolFromSmiles(prod_smi)
                fragments = apply_edits_to_mol(Chem.Mol(mol), edits=node.edit)

                tmp_frags = MultiElement(Chem.Mol(fragments)).mols
                fragments = Chem.Mol()
                for mol in tmp_frags:
                    fragments = Chem.CombineMols(fragments, mol)

                frag_graph = MultiElement(mol=Chem.Mol(fragments), rxn_class=rxn_class)

                directed = self.directed if hasattr(self, 'directed') else self.lg_directed

                frag_tensors, frag_scopes = pack_graph_feats(graph_batch=[frag_graph],
                                                             directed=directed,
                                                             return_graphs=False,
                                                             use_rxn_class=use_rxn_class)
                new_list.append(copy_node(node))
            except:
                continue
        return new_list

    def run_search(self, prod_smi, max_steps=6, rxn_class=None):
        with torch.no_grad():  # 无梯度计算
            node_list = self.run_edit_step(prod_smi, rxn_class=rxn_class)  # 预测的反应中心节点
            new_node_list = [copy_node(node) for node in node_list]  # 复制节点
            steps = 0

            # 每一个节点与产物作为输入，创建    特征组合？？？
            new_node_list = [self._create_lg_node(prod_smi, node, rxn_class=rxn_class) for node in new_node_list]

            while not check_nodes_complete(new_node_list) and steps <= max_steps:
                tmp_list = self.run_lg_step(prod_smi, new_node_list)
                new_node_list = [copy_node(node) for node in tmp_list]
                steps += 1

        new_node_list = self.keep_topk_nodes(new_node_list)
        return new_node_list

    def run_edit_step(self, prod_smi, rxn_class=None, **kwargs):
        node_list = self.get_topk_edits(prod_smi, rxn_class=rxn_class)
        return self.remove_invalid_nodes(prod_smi, node_list, rxn_class=rxn_class)

    def _create_lg_node(self, prod_smi, node, edit_2=None,rxn_class=None):  # 产物+断键信息集合node
        new_node = copy_node(node)
        mol = Chem.MolFromSmiles(prod_smi)

        use_rxn_class = False
        if rxn_class is not None:
            use_rxn_class = True

        # 设置有向图属性
        directed = self.directed if hasattr(self, 'directed') else self.lg_directed

        if not hasattr(self.model, 'encoder'):
            # print('encoder')   fff
            prod_graph = RxnElement(mol=Chem.Mol(mol), rxn_class=rxn_class)

            # 打包图的特征  prod_graph 的特征打包为张量
            prod_tensors, prod_scopes = pack_graph_feats([prod_graph],
                                                          directed=directed,
                                                          return_graphs=False,
                                                          use_rxn_class=use_rxn_class)
            prod_tensors = self.model.to_device(prod_tensors)

            # 编码产物向量
            prod_vecs, _ = self.model.lg_net.encoder(prod_tensors, prod_scopes)
            new_node.prod_vecs = prod_vecs.clone()  # new_node 的一个属性包含产物的特征


        fragments = apply_edits_to_mol(Chem.Mol(mol), new_node.edit)  # 得到合成子   【原来】

        # fragments = new_edits_split(5,6,15,16)

        # fragments_new1 = Chem.MolToSmiles(fragments)
        # fragments_new2 = []
        # fragments_new2.append(fragments_new1.split('.'))




        # print(fragments_new2[0][0])
        # new_edits_split(5,6,15,16)


        # 两个合成子的抽象信息  【x， y】
        tmp_frags = MultiElement(Chem.Mol(fragments)).mols

        fragments = Chem.Mol()
        for mol in tmp_frags:
            fragments = Chem.CombineMols(fragments, mol)   # 某一个合成子 + 整体合成子 ==> 该合成子的抽象表示

        frag_graph = MultiElement(mol=Chem.Mol(fragments), rxn_class=rxn_class)
        frag_tensors, frag_scopes = pack_graph_feats([frag_graph],
                                                      directed=directed,
                                                      return_graphs=False,
                                                      use_rxn_class=use_rxn_class)  #  特征融合
        frag_tensors = self.model.to_device(frag_tensors)

        if not hasattr(self.model, 'encoder'):
            frag_vecs, _ = self.model.lg_net.encoder(frag_tensors, frag_scopes)
        else:
            frag_vecs, _ = self.model.encoder(frag_tensors, frag_scopes)
        frag_vecs = torch.nn.utils.rnn.pad_sequence(frag_vecs, batch_first=True)
        assert len(frag_vecs.shape) == 3

        if not hasattr(self.model, 'config'):
            assert frag_vecs.shape[-1] == self.model.lg_net.config['mpn_size']
        else:
            assert frag_vecs.shape[-1] == self.model.config['mpn_size']
        new_node.num_fragments = frag_vecs.size(1)
        new_node.frag_vecs = frag_vecs.clone()
        return new_node ,fragments

    def add_lg_to_node(self, node):
        new_node = copy_node(node)
        if not new_node.node_complete:
            # print('new_node.frag_vecs[:, new_node.lg_idx]',new_node.frag_vecs[:, new_node.lg_idx])
            scores_lg, _ = self.model._compute_lg_step(graph_vecs=new_node.frag_vecs[:, new_node.lg_idx],
                                                    prod_vecs=new_node.prod_vecs.clone(),
                                                    prev_embed=new_node.prev_embed)

            scores_lg = F.log_softmax(scores_lg, dim=-1)  # 输出每个类别的对数概率


            if not hasattr(self.model, 'lg_vocab'):
                print('1', scores_lg.shape[-1])
                assert scores_lg.shape[-1] == len(self.model.lg_net.lg_vocab)
            else:
                print('2',scores_lg.shape[-1])
                assert scores_lg.shape[-1] == len(self.model.lg_vocab)

            #  topk_vals 包含前 k 个最高分数的张量; topk_idxs 包含对应于这些最高分数的索引的张量，这些索引可以用于后续的操作
            # print('scores_lg[0]', scores_lg[0])  # 保存了每一个词汇的分数，越靠近0越好
            topk_vals, topk_idxs = torch.topk(scores_lg[0], k=self.beam_width)
            # print('top',topk_idxs)


            # 把节点复制了5份
            new_list = [copy_node(new_node) for _ in range(self.beam_width)]
            # print(len(new_list))

            if hasattr(self.model, 'encoder'):
                for i_tensor, node in zip(*(topk_idxs, new_list)):  # 执行这里   topk_idxs 是前 k 个候选离去基团的索引（张量形式）   new_list 是当前步骤的候选节点列表
                    i = i_tensor.item() # 将张量形式的索引转换为 Python 整数
                    if isinstance(self.model.lg_embedding, nn.Linear):
                        node.prev_embed = self.model.lg_embedding(self.model.E_lg.index_select(index=i_tensor, dim=0))
                    else:
                        node.prev_embed = self.model.lg_embedding.index_select(index=i_tensor, dim=0)
                    node.add_lg(self.model.lg_vocab.get_elem(i), scores_lg[:, i].item())  # 索引离去基团以及对应的分数

            else:
                for i_tensor, node in zip(*(topk_idxs, new_list)):  #  使用 zip 将 topk_idxs 和 new_list 结合在一起，使得在每次迭代中可以同时访问最高分数的索引（i_tensor）和对应的节点（node）
                    i = i_tensor.item()
                    if isinstance(self.model.lg_net.lg_embedding, nn.Linear):
                        node.prev_embed = self.model.lg_net.lg_embedding(self.model.lg_net.E_lg.index_select(index=i_tensor, dim=0))
                    else:
                        node.prev_embed = self.model.lg_net.lg_embedding.index_select(index=i_tensor, dim=0)
                    node.add_lg(self.model.lg_net.lg_vocab.get_elem(i), scores_lg[:, i].item())
                    # a_sample = self.model.lg_net.lg_vocab.get_elem(i)
                    # a_score = scores_lg[:, i].item()

                    #  重要
                    # print('a_sample',a_sample)   # 预测的离去基团
                    # print('a_score',a_score)     # 其对应的分数

            assert all([node.frag_vecs is not None for node in new_list])
            return new_list

        else:
            if not hasattr(self.model, 'lg_vocab'):
                new_list = [copy_node(new_node) for _ in range(self.beam_width)]
            else:
                new_list = [copy_node(new_node) for _ in range(self.beam_width)]
            return new_list

    def run_lg_step(self, prod_smi, node_list):
        new_list = []
        assert all([node.frag_vecs is not None for node in node_list])
        for node in node_list:
            new_list.extend(self.add_lg_to_node(node))

        new_list = self.keep_topk_nodes(new_list)
        return new_list

    def get_edit_from_logits(self, mol, edit_logits, idx, val):
        if not hasattr(self.model, 'config'):
            config = self.model.edit_net.config
            toggles = self.model.edit_net.toggles
        else:
            config = self.model.config
            toggles = self.model.toggles

        if config['bs_outdim'] > 1:
            max_bond_idx = mol.GetNumBonds() * len(BOND_FLOATS)
        elif config['bs_outdim'] == 1:
            max_bond_idx = mol.GetNumBonds()

        if toggles.get('use_h_labels', False):

            if idx.item() < max_bond_idx:

                if config['bs_outdim'] > 1:
                    bond_logits = edit_logits[:mol.GetNumBonds() * len(BOND_FLOATS)]
                    bond_logits = bond_logits.reshape(mol.GetNumBonds(), len(BOND_FLOATS))
                    idx_tensor = torch.where(bond_logits == val)

                    idx_tensor = [indices[-1] for indices in idx_tensor]

                    bond_idx, bo_idx = idx_tensor[0].item(), idx_tensor[1].item()
                    a1 = mol.GetBondWithIdx(bond_idx).GetBeginAtom().GetAtomMapNum()
                    a2 = mol.GetBondWithIdx(bond_idx).GetEndAtom().GetAtomMapNum()

                    a1, a2 = sorted([a1, a2])
                    bo = mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble()
                    new_bo = BOND_FLOATS[bo_idx]

                    edit = f"{a1}:{a2}:{bo}:{new_bo}"

                elif config['bs_outdim'] == 1:
                    bond_idx = idx.item()
                    a1 = mol.GetBondWithIdx(bond_idx).GetBeginAtom().GetAtomMapNum()
                    a2 = mol.GetBondWithIdx(bond_idx).GetEndAtom().GetAtomMapNum()

                    a1, a2 = sorted([a1, a2])
                    bo = mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble()

                    edit = f"{a1}:{a2}:{bo}:{0.0}"

                else:
                    pass

            else:
                h_logits = edit_logits[max_bond_idx:]
                assert len(h_logits) == mol.GetNumAtoms()

                atom_idx = idx.item() - max_bond_idx
                a1 = mol.GetAtomWithIdx(atom_idx).GetAtomMapNum()

                edit = f"{a1}:{0}:{1.0}:{0.0}"

        else:

            if idx.item() == len(edit_logits) - 1:
                pass

            elif config['bs_outdim'] > 1:
                bond_logits = edit_logits[:mol.GetNumBonds() * len(BOND_FLOATS)]
                bond_logits = bond_logits.reshape(mol.GetNumBonds(), len(BOND_FLOATS))
                idx_tensor = torch.where(bond_logits == val)

                idx_tensor = [indices[-1] for indices in idx_tensor]

                bond_idx, bo_idx = idx_tensor[0].item(), idx_tensor[1].item()
                a1 = mol.GetBondWithIdx(bond_idx).GetBeginAtom().GetAtomMapNum()
                a2 = mol.GetBondWithIdx(bond_idx).GetEndAtom().GetAtomMapNum()

                a1, a2 = sorted([a1, a2])
                bo = mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble()
                new_bo = BOND_FLOATS[bo_idx]

                edit = f"{a1}:{a2}:{bo}:{new_bo}"

            elif config['bs_outdim'] == 1:
                bond_idx = idx.item()
                a1 = mol.GetBondWithIdx(bond_idx).GetBeginAtom().GetAtomMapNum()
                a2 = mol.GetBondWithIdx(bond_idx).GetEndAtom().GetAtomMapNum()

                a1, a2 = sorted([a1, a2])
                bo = mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble()

                edit = f"{a1}:{a2}:{bo}:{0.0}"

            else:
                pass

        return edit

class EditSearch(BeamSearch):

    def get_topk_edits(self, prod_smi, rxn_class=None, **kwargs):
        mol = Chem.MolFromSmiles(prod_smi)
        node_list = [BeamNode(mol=Chem.Mol(mol)) for _ in range(self.beam_width)]

        use_rxn_class = False
        if rxn_class is not None:
            use_rxn_class = True

        prod_graph = RxnElement(mol=Chem.Mol(mol), rxn_class=rxn_class)
        directed = self.directed if hasattr(self, 'directed') else self.edit_directed

        if self.max_edits == 1:
            prod_tensors, prod_scopes = pack_graph_feats([prod_graph],
                                                         directed=directed,
                                                         return_graphs=False,
                                                         use_rxn_class=use_rxn_class)

            prod_tensors = self.model.to_device(prod_tensors)
            bg_inputs = tensorize_bond_graphs([prod_graph], directed=directed,
                                               use_rxn_class=use_rxn_class)
            bg_tensors, bg_scope = bg_inputs
            bg_tensors = self.model.to_device(bg_tensors)
            bg_inputs = (bg_tensors, bg_scope)
            prod_vecs, edit_logits, _ = self.model._compute_edit_logits(prod_tensors, prod_scopes,
                                                                        ha=None, bg_inputs=bg_inputs)
            edit_logits = edit_logits[0]

            edit_logits = F.log_softmax(edit_logits, dim=-1)
            k = min(len(edit_logits), self.beam_width)
            topk_vals, topk_idxs = torch.topk(edit_logits, k=k)

            for beam_idx, (topk_idx, val) in enumerate(zip(*(topk_idxs, topk_vals))):
                edit = self.get_edit_from_logits(mol=Chem.Mol(mol),
                                                 edit_logits=edit_logits,
                                                 idx=topk_idx, val=val)
                if not isinstance(edit, list):
                    edit = [edit]
                node_list[beam_idx].add_edit(edit=edit, edit_prob=val.item())
                if hasattr(self.model, 'encoder'):
                    node_list[beam_idx].prod_vecs = prod_vecs.clone()

            return [copy_node(node) for node in node_list]
        else:
            raise ValueError("Greater than 1 sequence length not supported yet.")


class LGSearch(BeamSearch):
    def run_search(self, prod_smi, edits, rxn_class=None, max_steps=6):
        self.model.eval()
        with torch.no_grad():
            prod_mol = Chem.MolFromSmiles(prod_smi)
            prod_graph = RxnElement(mol=prod_mol, rxn_class=rxn_class)

            use_rxn_class = False
            if rxn_class is not None:
                use_rxn_class = True

            directed = self.directed if hasattr(self, 'directed') else self.lg_directed

            prod_tensors, prod_scopes = pack_graph_feats([prod_graph],
                                                          directed=directed,
                                                          return_graphs=False,
                                                          use_rxn_class=use_rxn_class)
            prod_tensors = self.model.to_device(prod_tensors)

            if hasattr(self.model, 'encoder'):
                prod_vecs, _ = self.model.encoder(prod_tensors, prod_scopes)
            else:
                prod_vecs, _ = self.model.lg_net.encoder(prod_tensors, prod_scopes)

            node_list = [BeamNode(mol=prod_mol)]
            for node in node_list:
                node.prod_vecs = prod_vecs.clone()
                node.add_edit(edits, 0.0)

            steps = 0
            new_node_list = [self._create_lg_node(prod_smi, node, rxn_class=rxn_class) for node in node_list]

            while not check_nodes_complete(new_node_list) and steps <= max_steps:
                tmp_list = self.run_lg_step(prod_smi, new_node_list)
                new_node_list = [copy_node(node) for node in tmp_list]
                steps += 1

            new_node_list = self.keep_topk_nodes(new_node_list)
            #print([node.lg_groups for node in new_node_list])

        return new_node_list
