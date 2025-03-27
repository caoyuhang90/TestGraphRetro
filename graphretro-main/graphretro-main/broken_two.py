import numpy as np
import pandas as pd
import torch
import os
import argparse
from tqdm import tqdm
from rdkit import RDLogger, Chem
from rdkit.Chem import Draw
import yaml

from seq_graph_retro.utils.parse import get_reaction_info, extract_leaving_groups
from seq_graph_retro.utils.chem import apply_edits_to_mol,get_mol
from seq_graph_retro.utils.edit_mol import canonicalize, generate_reac_set
from seq_graph_retro.models import EditLGSeparate
from seq_graph_retro.search import BeamSearch, BeamNode

from fastapi import FastAPI # fast api
from fastapi.middleware.cors import CORSMiddleware # 跨域问题
import uvicorn  # 直接启动APP
import base64  # 编码和解码 Base64
import io
from pydantic import BaseModel,Field
from fastapi.responses import JSONResponse

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 允许所有源；可以指定特定源，例如 ["https://example.com"]
    allow_credentials=True,  # 允许携带凭证，如 Cookies
    allow_methods=["*"],  # 允许所有 HTTP 方法，如 GET、POST、PUT、DELETE 等
    allow_headers=["*"],  # 允许所有请求头
) # 跨域处理


lg = RDLogger.logger()
lg.setLevel(4)

############################# 路径定义
try:
    ROOT_DIR = os.environ["SEQ_GRAPH_RETRO"]  # /home/wuhexing/GraphRetro
    DATA_DIR = os.path.join(ROOT_DIR, "datasets", "uspto-50k")
    EXP_DIR = os.path.join(ROOT_DIR, "models")

except KeyError:
    ROOT_DIR = "./"
    DATA_DIR = os.path.join(ROOT_DIR, "datasets", "uspto-50k")
    EXP_DIR = os.path.join(ROOT_DIR, "models")

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DEFAULT_TEST_FILE = f"{DATA_DIR}/canonicalized_test.csv"

print('DATA_DIR', DATA_DIR, 'EXP_DIR', EXP_DIR)




def canonicalize_prod(p):
    pcanon = canonicalize(p)
    pmol = Chem.MolFromSmiles(pcanon)
    [atom.SetAtomMapNum(atom.GetIdx() + 1) for atom in pmol.GetAtoms()]
    p = Chem.MolToSmiles(pmol)
    return p


def map_value_to_category(value,width_nub):
    return (value//(width_nub**2))


def check_nodes_complete(node_list):
    for node in node_list:
        if not node.node_complete:
            return False
    return True


# 【编辑】模型
def load_edits_model(args):
    edits_step = args.edits_step  # edits_step =  epoch_156
    if edits_step is None:
        edits_step = "best_model"
    if "run" in args.edits_exp:  # SingleEdit_10-02-2021--08-44-37
        # This addition because some of the new experiments were run using wandb
        edits_loaded = torch.load(os.path.join(args.exp_dir, "wandb", args.edits_exp, "files", edits_step + ".pt"),
                                  map_location=DEVICE)
        with open(f"{args.exp_dir}/wandb/{args.edits_exp}/files/config.yaml", "r") as f:
            tmp_loaded = yaml.load(f, Loader=yaml.FullLoader)
        model_name = tmp_loaded['model']['value']
    else:
        edits_loaded = torch.load(os.path.join(args.exp_dir, args.edits_exp,
                                               "checkpoints", edits_step + ".pt"),
                                  map_location=DEVICE)  # ./models/SingleEdit_10-02-2021--08-44-37/checkpoints/epoch_156.pt
        model_name = args.edits_exp.split("_")[0]  # model_name = SingleEdit
    return edits_loaded, model_name  # edits_loaded包含模型主体


# 【离去】模型
def load_lg_model(args):
    lg_step = args.lg_step  # step_101951
    if lg_step is None:
        lg_step = "best_model"
    if "run" in args.lg_exp:
        # This addition because some of the new experiments were run using wandb
        lg_loaded = torch.load(os.path.join(args.exp_dir, "wandb", args.lg_exp, "files", lg_step + ".pt"),
                               map_location=DEVICE)
        with open(f"{args.exp_dir}/wandb/{args.lg_exp}/files/config.yaml", "r") as f:
            tmp_loaded = yaml.load(f, Loader=yaml.FullLoader)
        model_name = tmp_loaded['model']['value']
    else:
        lg_loaded = torch.load(os.path.join(args.exp_dir, args.lg_exp,
                                            "checkpoints", lg_step + ".pt"),
                               map_location=DEVICE)  # models/LGIndEmbed_18-02-2021--12-23-26/checkpoints/step_101951.pt
        model_name = args.lg_exp.split("_")[0]  # model_name = LGIndEmbed
    return lg_loaded, model_name


def load_models(num):
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default=DATA_DIR, help="Data directory")
    parser.add_argument("--exp_dir", default=EXP_DIR, help="Experiments directory.")
    parser.add_argument("--test_file", default=DEFAULT_TEST_FILE, help="Test file.")
    parser.add_argument("--edits_exp", default="SingleEdit_10-02-2021--08-44-37",help="Name of edit prediction experiment.")
    parser.add_argument("--edits_step", default="epoch_156",help="Checkpoint for edit prediction experiment.")
    parser.add_argument("--lg_exp", default="LGIndEmbed_18-02-2021--12-23-26",help="Name of synthon completion experiment.")
    parser.add_argument("--lg_step", default="step_101951",help="Checkpoint for synthon completion experiment.")
    #beam_search 参数设置
    parser.add_argument("--beam_width", default=num, type=int, help="Beam width")
    parser.add_argument("--use_rxn_class", action='store_true', help="Whether to use reaction class.")
    parser.add_argument("--rxn_class_acc", action="store_true",help="Whether to print reaction class accuracy.")
    args, unknown = parser.parse_known_args()   #用一个unknown忽略不认识的参数
    test_df = pd.read_csv(args.test_file)    #读入测试文件\
    edits_loaded, edit_net_name = load_edits_model(args)  # 加载图编辑预测模型
    lg_loaded, lg_net_name = load_lg_model(args)  # 加载离去基团选择模型

    # 从加载的模型中提取配置数据
    edits_config = edits_loaded["saveables"]
    lg_config = lg_loaded['saveables']
    lg_toggles = lg_config['toggles']

    if 'tensor_file' in lg_config:  # 这里是false
        print('ok')
        if not os.path.isfile(lg_config['tensor_file']):
            if not lg_toggles.get("use_rxn_class", False):
                tensor_file = os.path.join(args.data_dir, "train/h_labels/without_rxn/lg_inputs.pt")
            else:
                tensor_file = os.path.join(args.data_dir, "train/h_labels/with_rxn/lg_inputs.pt")
            lg_config['tensor_file'] = tensor_file

    rm = EditLGSeparate(edits_config=edits_config, lg_config=lg_config, edit_net_name=edit_net_name,
                        lg_net_name=lg_net_name, device=DEVICE)  # 添加配置信息
    rm.load_state_dict(edits_loaded['state'], lg_loaded['state'])  # 加载模型信息  加载模型的权重和状态

    rm.to(DEVICE)  # 将模型移动到指定的设备（CPU 或 GPU）
    rm.eval()  # 设置模型为评估模式
    n_matched = np.zeros(args.beam_width)  # beam_size,默认为10
    width_nub = args.beam_width
    print('束搜索宽度', args.beam_width)

    beam_model = BeamSearch(model=rm, beam_width=args.beam_width, max_edits=1)  # 加载束搜索模型

    return width_nub, beam_model



def two_breaken_new(smiles,a1,a2,a11,a22,width_nub,beam_model,b1=1.0, b2=0,b11=1.0,b22=0,first_num=5,second_num=8):
    print('执行：broken_tow_bonds')
    rxn = []
    rxn1 = []
    rxn2 = []
    for id, smile in enumerate(smiles):
        p = canonicalize_prod(smile)  # 清理 + 重新编码
        rxn.append(p)  # 将产物加入到列表中
        mol = Chem.MolFromSmiles(p)
        # node_list=beam_model.run_edit_step(p)  # 预测反应中心
        node_list = []
        node_list.append(BeamNode(mol=Chem.Mol(mol)))  # [BeamNode(mol=Chem.Mol(mol))
        node_list[0].edit = [f'{a1}:{a2}:{b1}:{b2}']
        edit_2 = f'{a11}:{a22}:{b11}:{b22}'
        # print('模型预测的断键：',node_list[0].edit)
        new_node_list = []
        step = 0
        new_node, fragments_new2 = beam_model._create_lg_node(p, node_list[0],
                                                              rxn_class=None)  # 增加了两个特征向量（一个关于合成子frag_vecs，一个关于产物prod_vecs），合成子个数
        new_node_list.append(new_node)

        num_fragments = new_node_list[0].num_fragments  # 合成子个数

        # 第二次拆分的 【键】
        a1, a2, b1, b2 = edit_2.split(':')
        a1, a2, b1, b2 = int(a1), int(a2), float(b1), float(b2)

        while not check_nodes_complete(new_node_list) and step <= 6:
            tmp_list = []
            assert all([node.frag_vecs is not None for node in new_node_list])
            for node in new_node_list:
                tmp_list.extend(beam_model.add_lg_to_node(node))  # 预测离去基团，并添加至node中
            new_node_list = tmp_list
            step += 1

        step = 0
        tmp_list1 = tmp_list[:25]
        tmp_list1 = beam_model.keep_topk_nodes(tmp_list1)  # 选最靠谱的 5 个
        tmp_list2 = []  # 从第一次合成子

        tem_list1_scores = [node.prob for node in tmp_list1]  # 第一次预测的分数

        frag_1 = []
        frag_2 = []
        frag_FirstToSecond = []  # 保存第二次要拆的合成子的smiles，带编号

        for beam_idx, node in enumerate(tmp_list1):
            pred_edit = node.edit
            pred_label = node.lg_groups
            print(beam_idx, "预测的离去基团的结果", pred_label)

            if isinstance(pred_edit, list):  # 检查pre_edit是否为一个列表，如果是，执行下面语句
                pred_edit = pred_edit[0]
            try:
                pred_set, reac_smi_no_canonicalize = generate_reac_set(p, pred_edit, pred_label,
                                                                       verbose=False)  # 生成预测的反应物集合
            except BaseException as e:
                print(e, flush=True)
                pred_set = None

            pred_list = list(pred_set)  # 反应集转换为列表，并打印结果
            print('pred_list', pred_list)

            # 判断第二次断键所在的合成子：读取合成子们的原子标签，然后和指定标签对比
            frag_1 = reac_smi_no_canonicalize.split('.')  # frag_1保存合成子的smiles
            frag_2 = [Chem.MolFromSmiles(smi) for smi in frag_1]
            try:  # 生成的结果转化为mol后为None，跳过
                frag_amaps = [set([atom.GetAtomMapNum() for atom in mol.GetAtoms()])
                              for mol in frag_2]
            except BaseException as e:
                print(e, flush=True)
                continue
            if new_node_list[0].num_fragments == 1:
                frag_FirstToSecond.append(reac_smi_no_canonicalize)
                mol = Chem.MolFromSmiles(reac_smi_no_canonicalize)
                node_list.append(BeamNode(mol=Chem.Mol(mol)))
                node_list[1].edit = [edit_2]
                new_node, _ = beam_model._create_lg_node(reac_smi_no_canonicalize, node_list[1], rxn_class=None)
                tmp_list2.append(new_node)
            elif new_node_list[0].num_fragments == 2:
                if a1 in frag_amaps[0] and a2 in frag_amaps[0]:
                    print('frag_amaps[0]')
                    sml = frag_1[0]
                    frag_FirstToSecond.append(sml)
                    mol = Chem.MolFromSmiles(sml)
                    node_list.append(BeamNode(mol=Chem.Mol(mol)))
                    node_list[1].edit = [edit_2]
                    new_node, _ = beam_model._create_lg_node(sml, node_list[1], rxn_class=None)
                    tmp_list2.append(new_node)
                elif a1 in frag_amaps[1] and a2 in frag_amaps[1]:
                    print('frag_amaps[1]')
                    # sml = canonicalize_prod(fragments_new2[0][1])
                    sml = frag_1[1]
                    frag_FirstToSecond.append(sml)
                    mol = Chem.MolFromSmiles(sml)
                    node_list.append(BeamNode(mol=Chem.Mol(mol)))
                    node_list[1].edit = [edit_2]
                    new_node, _ = beam_model._create_lg_node(sml, node_list[1], rxn_class=None)
                    tmp_list2.append(new_node)
                # print('第二次断键选择的第一次的合成子sml', sml)

            if len(pred_list) > 1:
                rxn1.append(pred_list[0] + '.' + pred_list[1])
            else:
                rxn1.append(pred_list[0])

        # frag_FirstToSecond 存放对应前体
        frag_FirstToSecond_new = [canonicalize(x) for x in frag_FirstToSecond]

        if new_node_list[0].num_fragments == 1:
            rxn1 = [None]
        elif new_node_list[0].num_fragments == 2:
            for id, x in enumerate(frag_FirstToSecond):
                x = canonicalize(x)
                y = rxn1[id].split('.')
                y.remove(x)
                rxn1[id] = y[0]


        # 束搜索匹配离去基团，关键参数 tmp_list2    开始有5个，然后扩展为25个，再扩展为125个
        while not check_nodes_complete(tmp_list2) and step <= 6:
            tmp_list = []
            assert all([node.frag_vecs is not None for node in tmp_list2])
            for index, node in enumerate(tmp_list2):
                tmp_list.extend(beam_model.add_lg_to_node(node))  # 重要 预测离去基团，并添加至node中
            tmp_list2 = tmp_list
            step += 1

        tmp_list2 = beam_model.keep_topk_nodes_new(tmp_list2, second_num)  # 选最靠谱的 second_num 个   现在的tmp_list2元组表示

        tmp_list2 = [[map_value_to_category(node[0], width_nub), node[1]] for node in tmp_list2]

        tmp_list2_scores = [(node[0], node[1].prob) for node in tmp_list2]
        # 得分处理  第一次预测得分与第二次得分的综合得分
        CompositeScores = []  # 存放分数
        scores_list = []
        scores_list_true = []
        for id, node in enumerate(tmp_list2_scores):
            CompositeScores.append((node[0], node[1] * tem_list1_scores[node[0]]))
        # 分数重置tmp_list2
        for id, node in enumerate(tmp_list2):
            composite_score = CompositeScores[id][1]
            scores_list.append(composite_score)
            node[1].prob = composite_score
        # 再排序
        sorted_list_2 = sorted(tmp_list2, key=lambda x: x[1].prob, reverse=False)

        for id, node in enumerate(sorted_list_2):
            node.append(frag_FirstToSecond[node[0]])
            # node.append(rxn1[node[0]])
            if  rxn1[0] is not None:
                node.append(rxn1[node[0]])
            else:
                node.append(None)



    return sorted_list_2,num_fragments

# 分子式子，键1，键2，第一次束宽，第二次束宽
def two_breaken_new001(sml,edit_1,edit_2,beam_num01=5,beam_num02=5):

    beam_num01 = beam_num01
    width_nub, beam_model = load_models(beam_num01)  # 模型加载

    beam_num02 = beam_num02
    smile = [sml]
    edit_1 = edit_1
    edit_2 = edit_2

    a1, a2 = edit_1.split(":")
    a1, a2 = int(a1), int(a2)
    a11, a22 = edit_2.split(":")
    a11, a22 = int(a11), int(a22)
    list_001 = []
    list_002 = []

    rxn = []
    rxn1 = []
    rxn2 = []

    list_001, num_fragments = two_breaken_new(smile, a1, a2, a11, a22, width_nub=width_nub, beam_model=beam_model,
                                              second_num=beam_num02)
    list_002, _ = two_breaken_new(smile, a11, a22, a1, a2, width_nub=width_nub, beam_model=beam_model,
                                  second_num=beam_num02)

    list_000 = list_001 + list_002

    sorted_list_2 = sorted(list_000, key=lambda x: x[1].prob, reverse=False)[:beam_num02]

    for beam_idx, node in enumerate(sorted_list_2):  # tmp_list2 是元组
        pred_edit = node[1].edit
        pred_label = node[1].lg_groups

        sml = node[2]

        if isinstance(pred_edit, list):  # 检查pre_edit是否为一个列表，如果是，执行下面语句
            pred_edit = pred_edit[0]

        print(beam_idx, "离去基团的标签预测结果", pred_label)

        try:
            # rxn2.append()
            pred_set, _ = generate_reac_set(sml, pred_edit, pred_label, verbose=False)  # 重要  生成预测的反应集
        except BaseException as e:
            print(e, flush=True)
            pred_set = None

        # sml = Chem.MolToSmiles(sml)
        print('sml', sml)
        print(beam_idx, "离去基团的标签预测结果", pred_label)
        pred_list = list(pred_set)  # 反应集转换为列表，并打印结果

        print('pred_list', pred_list)
        if len(pred_list) > 1:
            rxn2.append(pred_list[0] + '.' + pred_list[1])
        else:
            rxn2.append(pred_list[0])

        if num_fragments == 1:
            rxn.append(rxn2[beam_idx])
        elif num_fragments == 2:
            rxn.append(node[3] + '.' + rxn2[beam_idx])

    print('rxn', rxn)
    mol_rxn = [Chem.MolFromSmiles(sml) for sml in rxn]
    img = Draw.MolsToGridImage(mol_rxn, molsPerRow=3, subImgSize=(800, 800))
    img.save("aab_duanjian/000结果.png")




# smile = 'CCNC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1/C=C/[C@H](CCC2=CC=CC=C2)O)O)O'
# edit_1 = '9:10'
# edit_2 = '14:15'

smile = 'O=C1[C@@H](C2C=CC=CC=2)O[C@@H](C(C)(C)C)O1'
edit_1 = '5:16'
edit_2 = '6:7'

two_breaken_new001(smile,edit_1,edit_2,5,10)

# beam_num01 = 5
# beam_num02 = 10
# width_nub, beam_model = load_models(beam_num01)
# smile = ['CCNC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1/C=C/[C@H](CCC2=CC=CC=C2)O)O)O']
# edit_1 = '9:10'
# edit_2 = '14:15'
#
# # smile = ['O=C1[C@@H](C2C=CC=CC=2)O[C@@H](C(C)(C)C)O1']
# # edit_1 = '5:16'
# # edit_2 = '6:7'
#
# a1, a2 = edit_1.split(":")
# a1, a2 = int(a1), int(a2)
# a11, a22 = edit_2.split(":")
# a11, a22 = int(a11), int(a22)
# list_001 = []
# list_002 = []
#
# rxn=[]
# rxn1=[]
# rxn2=[]
#
# list_001 ,num_fragments = two_breaken_new(smile,a1,a2,a11,a22,width_nub=width_nub,beam_model=beam_model,second_num = beam_num02)
# list_002 , _ = two_breaken_new(smile,a11,a22,a1,a2,width_nub=width_nub,beam_model=beam_model,second_num = beam_num02)
#
# list_000 = list_001+list_002
#
# sorted_list_2 = sorted(list_000, key=lambda x: x[1].prob, reverse=False)[:beam_num02]
#
# for beam_idx, node in enumerate(sorted_list_2):  # tmp_list2 是元组
#     pred_edit = node[1].edit
#     pred_label = node[1].lg_groups
#
#     sml = node[2]
#
#     if isinstance(pred_edit, list):  # 检查pre_edit是否为一个列表，如果是，执行下面语句
#         pred_edit = pred_edit[0]
#
#     print(beam_idx, "离去基团的标签预测结果", pred_label)
#
#     try:
#         # rxn2.append()
#         pred_set, _ = generate_reac_set(sml, pred_edit, pred_label, verbose=False)  # 重要  生成预测的反应集
#     except BaseException as e:
#         print(e, flush=True)
#         pred_set = None
#
#     # sml = Chem.MolToSmiles(sml)
#     print('sml', sml)
#     print(beam_idx, "离去基团的标签预测结果", pred_label)
#     pred_list = list(pred_set)  # 反应集转换为列表，并打印结果
#
#     print('pred_list', pred_list)
#     if len(pred_list) > 1:
#         rxn2.append(pred_list[0] + '.' + pred_list[1])
#     else:
#         rxn2.append(pred_list[0])
#
#     if num_fragments == 1:
#         rxn.append(rxn2[beam_idx])
#     elif num_fragments == 2:
#         rxn.append(node[3] + '.' + rxn2[beam_idx])
#
# print('rxn',rxn)
# mol_rxn = [Chem.MolFromSmiles(sml) for sml in rxn]
# img=Draw.MolsToGridImage(mol_rxn, molsPerRow=3, subImgSize=(800, 800))
# img.save("aab_duanjian/000结果.png")





print('ok')
