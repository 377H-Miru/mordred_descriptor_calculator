import sys, os
import pytest
import pandas as pd
import numpy as np
from rdkit import Chem

# 本体ファイルのパスを明示的に追加
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    from calculate_descriptors import setup_mordred_calculator, calc_conjugation_features
except ImportError as e:
    # デバッグ用にパス情報を出力
    print(f"DEBUG: sys.path = {sys.path}")
    print(f"DEBUG: Current dir = {os.getcwd()}")
    print(f"DEBUG: Files in root = {os.listdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))}")
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
