from rdkit import Chem

smiles = "[CH3:1][O:2][c:3]1[cH:4][c:5]2[c:6]([cH:7][c:8]1[O:9][CH3:10])-[c:11]1[cH:12][c:13]([Cl:1000])[nH:23][c:24](=[O:25])[n:26]1[CH2:27][CH2:28]2"

try:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Failed to parse SMILES, returned None.")
    else:
        print("Successfully parsed SMILES.")
except Exception as e:
    print(f"Error processing SMILES: {smiles}, Error: {e}")