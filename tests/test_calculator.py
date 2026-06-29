import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred_descriptor_calculator.core import setup_mordred_calculator, preprocess_molecule
from mordred_descriptor_calculator.conjugation import calc_conjugation_features

def test_mordred_calculator_setup():
    calc = setup_mordred_calculator(ignore_3d=True)
    assert calc is not None
    assert len(calc.descriptors) > 0

def test_preprocess_molecule_2d_vs_3d():
    mol_2d, can_2d, err_2d = preprocess_molecule("CCO", include_3d=False)
    assert mol_2d is not None
    assert mol_2d.GetNumConformers() == 0
    assert err_2d is None

    mol_3d, can_3d, err_3d = preprocess_molecule("CCO", include_3d=True, optimize=False)
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() > 0
    assert err_3d is None

def test_conjugation_benzene_and_ethane():
    mol_b = Chem.MolFromSmiles("c1ccccc1")
    res_b = calc_conjugation_features(mol_b)
    assert res_b["Conjugation_Count"] >= 1
    assert res_b["Conjugation_MaxAtomCount"] == 6
    assert res_b["Conjugation_Error"] == ""

    mol_e = Chem.MolFromSmiles("CC")
    res_e = calc_conjugation_features(mol_e)
    assert res_e["Conjugation_Count"] == 0
    assert res_e["Conjugation_MaxAtomCount"] == 0
    assert res_e["Conjugation_Error"] == ""

def test_conjugation_failure_distinction(monkeypatch):
    from mordred_descriptor_calculator import conjugation
    def mock_connected_components(graph):
        raise RuntimeError("Simulated Graph Error")

    monkeypatch.setattr(conjugation.nx, "connected_components", mock_connected_components)
    mol = Chem.MolFromSmiles("c1ccccc1")
    res = calc_conjugation_features(mol)
    assert res["Conjugation_Error"] == "Simulated Graph Error"

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
