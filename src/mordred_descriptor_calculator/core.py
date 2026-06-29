import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .compat import patch_numpy_for_mordred_compat
patch_numpy_for_mordred_compat()

from mordred import Calculator, descriptors
from mordred.PathCount import PathCount

def preprocess_molecule(smiles, seed=42, optimize=True, standardize=True, desalt=True, include_3d=False):
    """
    Preprocess a single SMILES string through parse, standardize, desalt, add_hydrogen, embed_3d, optimize_3d stages.
    Returns: (mol, canonical_smiles, error_info)
    """
    if pd.isna(smiles) or str(smiles).strip() == "":
        return None, "", {'stage': 'parse', 'error_type': 'EmptySMILES', 'error_message': 'SMILES string is empty or NaN'}
    
    clean_smiles = str(smiles).strip()
    
    # 1. Parse stage
    try:
        mol = Chem.MolFromSmiles(clean_smiles)
        if mol is None:
            return None, "", {'stage': 'parse', 'error_type': 'ParseError', 'error_message': 'RDKit unable to parse SMILES'}
    except Exception as e:
        return None, "", {'stage': 'parse', 'error_type': 'ParseException', 'error_message': str(e)}

    # 2. Standardize stage
    if standardize:
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            mol = rdMolStandardize.Cleanup(mol)
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)
        except Exception as e:
            return None, "", {'stage': 'standardize', 'error_type': 'StandardizeError', 'error_message': str(e)}

    # 3. Desalt stage
    if desalt:
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            lfc = rdMolStandardize.LargestFragmentChooser()
            mol = lfc.choose(mol)
        except Exception as e:
            try:
                frags = Chem.GetMolFrags(mol, asMols=True)
                if len(frags) > 1:
                    mol = max(frags, key=lambda m: m.GetNumAtoms())
            except Exception as ex:
                return None, "", {'stage': 'desalt', 'error_type': 'DesaltError', 'error_message': str(ex)}

    canonical_smiles = Chem.MolToSmiles(mol)

    if not include_3d:
        return mol, canonical_smiles, None

    # 4. Add Hydrogen stage
    try:
        mol = Chem.AddHs(mol)
    except Exception as e:
        return None, canonical_smiles, {'stage': 'add_hydrogen', 'error_type': 'AddHydrogenError', 'error_message': str(e)}

    # 5. Embed 3D stage
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        res = AllChem.EmbedMolecule(mol, params)
        if res == -1:
            return None, canonical_smiles, {'stage': 'embed_3d', 'error_type': 'Embed3DError', 'error_message': '3D coordinate embedding failed'}
    except Exception as e:
        return None, canonical_smiles, {'stage': 'embed_3d', 'error_type': 'Embed3DException', 'error_message': str(e)}

    # 6. Optimize 3D stage
    if optimize:
        try:
            opt_res = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            if opt_res == -1:
                opt_uff = AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                if opt_uff == -1:
                    return None, canonical_smiles, {'stage': 'optimize_3d', 'error_type': 'Optimize3DError', 'error_message': '3D structure optimization (MMFF & UFF) failed'}
        except Exception as e:
            return None, canonical_smiles, {'stage': 'optimize_3d', 'error_type': 'Optimize3DException', 'error_message': str(e)}

    return mol, canonical_smiles, None

def preprocess_worker(args):
    """
    Worker tuple wrapper for parallel execution.
    args: (smiles, row_props, seed, optimize, standardize, desalt, include_3d, row_idx, mol_id)
    """
    smiles, props, seed, optimize, standardize, desalt, include_3d, row_idx, mol_id = args
    mol, canonical_smiles, err = preprocess_molecule(
        smiles=smiles, seed=seed, optimize=optimize,
        standardize=standardize, desalt=desalt, include_3d=include_3d
    )
    if mol is not None:
        return mol, canonical_smiles, props, None
    else:
        err_entry = {
            'row_index': row_idx,
            'ID': mol_id,
            'input_smiles': str(smiles) if smiles is not None else "",
            'stage': err.get('stage', 'unknown'),
            'error_type': err.get('error_type', 'Error'),
            'error_message': err.get('error_message', '')
        }
        return None, canonical_smiles, props, err_entry

def setup_mordred_calculator(ignore_3d=True, descriptor_set="all"):
    """
    Build and configure Mordred Calculator.
    """
    if descriptor_set == "2d":
        calc = Calculator(descriptors, ignore_3D=True)
    elif descriptor_set == "3d":
        calc = Calculator(descriptors, ignore_3D=False)
    else:
        calc = Calculator(descriptors.all, ignore_3D=ignore_3d)
        
    for i in range(1, 51):
        try:
            calc.register(PathCount(order=i, pi=False))
            calc.register(PathCount(order=i, pi=True))
        except Exception:
            pass
    return calc
