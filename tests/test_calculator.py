import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred_descriptor_calculator.core import setup_mordred_calculator, preprocess_molecule
from mordred_descriptor_calculator.conjugation import calc_conjugation_features

def test_mordred_calculator_setup():
    calc = setup_mordred_calculator(ignore_3d=True, descriptor_set="2d")
    assert calc is not None
    assert len(calc.descriptors) > 0

def test_conjugation_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    res = calc_conjugation_features(mol)
    assert res["Conjugation_Count"] >= 1
    assert res["Conjugation_MaxAtomCount"] == 6

def test_conjugation_ethane():
    mol = Chem.MolFromSmiles("CC")
    res = calc_conjugation_features(mol)
    assert res["Conjugation_Count"] == 0
    assert res["Conjugation_MaxAtomCount"] == 0

def test_invalid_smiles_parse_error():
    mol, can_smiles, err = preprocess_molecule("INVALID_SMILES")
    assert mol is None
    assert err is not None
    assert err["stage"] == "parse"
    assert err["error_type"] == "ParseError"

def test_embed_3d_failure_logged(monkeypatch):
    def mock_embed_molecule(mol, params):
        return -1

    monkeypatch.setattr(AllChem, "EmbedMolecule", mock_embed_molecule)
    
    mol, can_smiles, err = preprocess_molecule("CCO", include_3d=True, optimize=False)
    assert mol is None
    assert err is not None
    assert err["stage"] == "embed_3d"
    assert err["error_type"] == "Embed3DError"

def test_optimize_3d_failure_logged(monkeypatch):
    def mock_mmff_optimize(mol, maxIters):
        return -1
    def mock_uff_optimize(mol, maxIters):
        return -1

    monkeypatch.setattr(AllChem, "MMFFOptimizeMolecule", mock_mmff_optimize)
    monkeypatch.setattr(AllChem, "UFFOptimizeMolecule", mock_uff_optimize)

    mol, can_smiles, err = preprocess_molecule("CCO", include_3d=True, optimize=True)
    assert mol is None
    assert err is not None
    assert err["stage"] == "optimize_3d"
    assert err["error_type"] == "Optimize3DError"
