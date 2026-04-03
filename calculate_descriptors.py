import os
import sys
import warnings
import argparse
import json
import concurrent.futures
import numpy as np

# Monkey-patch for mordred compatibility with modern numpy (CRITICAL)
if not hasattr(np, 'product'): np.product = np.prod

import networkx as nx
import pandas as pd
from mordred import Calculator, descriptors
from mordred.error import DuplicatedDescriptorName
from mordred.PathCount import PathCount
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

def preprocess_worker(args):
    """
    Worker function for parallel preprocessing.
    Args: (smiles, seed, optimize)
    Returns: (clean_mol, props_dict, error_info)
    """
    smiles, props, seed, optimize = args
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None: return None, None, {"smiles": smiles, "reason": "RDKit Parse Error"}
        
        # Standardization
        mol = rdMolStandardize.Cleanup(mol)
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        
        # 3D Conformer Generation with DETERMINISTIC SEED
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed # Fixed seed per molecule
        if AllChem.EmbedMolecule(mol, params) == -1:
            return None, None, {"smiles": smiles, "reason": "3D Embedding Failed"}
            
        if optimize:
            if AllChem.MMFFOptimizeMolecule(mol, maxIters=200) == -1:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                
        return mol, props, None
    except Exception as e:
        return None, None, {"smiles": smiles, "reason": f"Exception: {str(e)}"}

def calc_conjugation_features(mol):
    """Custom features for pi-conjugated systems."""
    res = {"Conjugation_Count": 0, "Conjugation_MaxAtomCount": 0, "Conjugation_MaxLength": 0, "Conjugation_BLA": np.nan, "Conjugation_GraphEnergy": 0.0}
    if mol is None or mol.GetNumConformers() == 0: return res
    try:
        conjugated_bonds = [b for b in mol.GetBonds() if b.GetIsConjugated()]
        if not conjugated_bonds: return res
        conf = mol.GetConformer()
        graph = nx.Graph()
        for b in conjugated_bonds:
            u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            w = conf.GetAtomPosition(u).Distance(conf.GetAtomPosition(v))
            graph.add_edge(u, v, weight=w)
        systems = list(nx.connected_components(graph))
        if not systems: return res
        largest = max(systems, key=len)
        sub = graph.subgraph(largest)
        res["Conjugation_Count"] = len(systems)
        res["Conjugation_MaxAtomCount"] = len(largest)
        if sub.number_of_nodes() > 1:
            try: res["Conjugation_MaxLength"] = nx.diameter(sub)
            except: pass
        lengths = [d['weight'] for u, v, d in sub.edges(data=True)]
        if lengths: res["Conjugation_BLA"] = np.std(lengths)
        adj = nx.to_numpy_array(sub, weight=None)
        res["Conjugation_GraphEnergy"] = np.sum(np.abs(np.linalg.eigvalsh(adj)))
        return res
    except: return res

def setup_mordred_calculator():
    calc = Calculator(descriptors.all, ignore_3D=False)
    for i in range(1, 51):
        try:
            calc.register(PathCount(order=i, pi=False))
            calc.register(PathCount(order=i, pi=True))
        except: pass
    return calc

def main():
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser(description="Deterministic & Reproducible Descriptor Calculator.")
    parser.add_argument('--config', help='JSON configuration file.')
    parser.add_argument('--seed', type=int, default=42, help='Base random seed.')
    parser.add_argument('--no-optimize', action='store_true', help='Skip energy minimization.')
    args = parser.parse_args()

    if not args.config: sys.exit("Error: Please provide --config file.")
    with open(args.config, 'r') as f: config = json.load(f)
    in_path, out_path = config['input_path'], config['output_path']
    smiles_col = config.get('smiles_col', 'smiles')

    df_in = pd.read_csv(in_path)
    # Prepare deterministic tasks (seed is incremented per row to ensure variety yet reproducibility)
    tasks = [(row[smiles_col], row.to_dict(), args.seed + i, not args.no_optimize) 
             for i, row in df_in.iterrows()]

    print(f"Parallel preprocessing {len(tasks)} compounds with seed {args.seed}...")
    
    mols, props, errors = [], [], []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(preprocess_worker, tasks), total=len(tasks), desc="Preprocessing"))

    for mol, p, err in results:
        if mol:
            mols.append(mol)
            props.append(p)
        else:
            errors.append(err)

    if mols:
        calc = setup_mordred_calculator()
        print("\nStarting Mordred descriptor calculation...")
        m_df = calc.pandas(mols, nproc=os.cpu_count(), quiet=False)
        m_df = m_df.apply(pd.to_numeric, errors='coerce')
        
        with concurrent.futures.ProcessPoolExecutor() as exec:
            c_list = list(exec.map(calc_conjugation_features, mols))
        c_df = pd.DataFrame(c_list)
        
        final_df = pd.concat([pd.DataFrame(props), m_df, c_df], axis=1)
        final_df.to_csv(out_path, index=False)
        print(f"Results saved to {out_path}")
    
    if errors:
        pd.DataFrame(errors).to_csv("errors.log", index=False)
        print(f"Warning: {len(errors)} molecules failed. See errors.log")

if __name__ == "__main__":
    main()
