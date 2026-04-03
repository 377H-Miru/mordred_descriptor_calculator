import numpy as np
if not hasattr(np, "product"): np.product = np.prod
import sys, os
import pytest
from rdkit import Chem

# 1. ワークスペースの構造を徹底的にデバッグ表示
print(f"\n--- DEBUG START ---")
print(f"Current Working Directory: {os.getcwd()}")
print(f"Script File: {__file__}")
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
print(f"Parent Directory: {parent_dir}")
print(f"Files in Parent: {os.listdir(parent_dir)}")

# 2. 親ディレクトリをパスの最優先に追加
sys.path.insert(0, parent_dir)

# 3. インポート試行
try:
    import calculate_descriptors
    print("SUCCESS: calculate_descriptors imported.")
    # 関数が存在するか確認
    if hasattr(calculate_descriptors, 'setup_mordred_calculator'):
        print("SUCCESS: setup_mordred_calculator found.")
        from calculate_descriptors import setup_mordred_calculator, calc_conjugation_features
    else:
        print(f"FAILURE: setup_mordred_calculator NOT found in {calculate_descriptors.__file__}")
        # ファイルの中身を数行表示して確認
        with open(calculate_descriptors.__file__, 'r') as f:
            print("--- FILE PEEK ---")
            print("".join(f.readlines()[:20]))
        raise ImportError("Function missing in the imported module.")
except ImportError as e:
    print(f"ERROR: Import failed. {e}")
    raise e

def test_mordred_setup():
    calc = setup_mordred_calculator()
    assert len(calc) > 1000

def test_custom_conjugation_descriptors():
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=42)
    results = calc_conjugation_features(mol)
    assert results["Conjugation_Count"] == 1
    assert results["Conjugation_MaxAtomCount"] == 6
    assert results["Conjugation_GraphEnergy"] > 0

def test_handling_non_conjugated_molecule():
    mol = Chem.MolFromSmiles("CC")
    results = calc_conjugation_features(mol)
    assert results["Conjugation_Count"] == 0
    assert results["Conjugation_MaxAtomCount"] == 0
    assert results["Conjugation_GraphEnergy"] == 0.0

def test_batch_processing_logic():
    calc = setup_mordred_calculator()
    mol = Chem.MolFromSmiles("CCO")
    mordred_res = calc(mol).asdict()
    assert len(mordred_res) > 0
    any_value = any(not isinstance(v, Exception) for v in mordred_res.values())
    assert any_value
