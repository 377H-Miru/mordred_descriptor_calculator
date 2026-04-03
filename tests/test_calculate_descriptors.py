import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import pytest
import pandas as pd
import numpy as np
from rdkit import Chem
from calculate_descriptors import setup_mordred_calculator, calc_conjugation_features

def test_mordred_setup():
    """Mordred計算機が正しく初期化されるか検証"""
    calc = setup_mordred_calculator()
    assert len(calc) > 1000  # 大量の記述子が登録されているはず

def test_custom_conjugation_descriptors():
    """独自実装の共役系記述子の計算ロジックを検証"""
    # ベンゼン（共役系1つ、原子6個）
    mol = Chem.MolFromSmiles("c1ccccc1")
    # 3D座標が必要なため生成
    mol = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=42)
    
    results = calc_conjugation_features(mol)
    
    assert results["Conjugation_Count"] == 1
    assert results["Conjugation_MaxAtomCount"] == 6
    assert results["Conjugation_GraphEnergy"] > 0

def test_handling_non_conjugated_molecule():
    """共役系を持たない分子（例：エタン）での動作を検証"""
    mol = Chem.MolFromSmiles("CC")
    results = calc_conjugation_features(mol)
    
    assert results["Conjugation_Count"] == 0
    assert results["Conjugation_MaxAtomCount"] == 0
    assert results["Conjugation_GraphEnergy"] == 0.0

def test_batch_processing_logic():
    """コアとなる計算関数の正常系のみを検証"""
    calc = setup_mordred_calculator()
    mol = Chem.MolFromSmiles("CCO")
    # Mordredは2Dでも動作する記述子が多い
    mordred_res = calc(mol)
    
    assert not mordred_res.fill_missing().isnull().any()
