import os
import sys
import warnings
import argparse
import json
import concurrent.futures
import numpy as np
import networkx as nx
import pandas as pd
from mordred import Calculator, descriptors
from mordred.error import DuplicatedDescriptorName
from mordred.PathCount import PathCount
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

# Monkey-patch for numpy 1.25+ compatibility
if not hasattr(np, 'product'): np.product = np.prod

def preprocess_molecule(mol, seed=42):
    """Standardize molecule and generate 3D conformer with fixed seed."""
    if mol is None: return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed # FIXED SEED FOR REPRODUCIBILITY
        if AllChem.EmbedMolecule(mol, params) == -1:
            return None # Skip if 3D embedding fails
        return mol
    except: return None

def calc_conjugation_features(mol):
    """Custom features for pi-conjugated systems with error handling."""
    res = {"Conjugation_Count": 0, "Conjugation_MaxAtomCount": 0, "Conjugation_MaxLength": 0, "Conjugation_BLA": np.nan, "Conjugation_GraphEnergy": 0.0}
    if mol is None: return res
    try:
        conjugated_bonds = [b for b in mol.GetBonds() if b.GetIsConjugated()]
        if not conjugated_bonds: return res
        
        has_3d = mol.GetNumConformers() > 0
        conf = mol.GetConformer() if has_3d else None
        
        graph = nx.Graph()
        for b in conjugated_bonds:
            u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            w = conf.GetAtomPosition(u).Distance(conf.GetAtomPosition(v)) if has_3d else 1.0
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
        if has_3d and sub.number_of_edges() > 0:
            res["Conjugation_BLA"] = np.std(list(nx.get_edge_attributes(sub, 'weight').values()))
            
        adj = nx.to_numpy_array(sub, weight=None)
        res["Conjugation_GraphEnergy"] = np.sum(np.abs(np.linalg.eigvalsh(adj)))
        return res
    except: return res

def setup_mordred_calculator():
    """Create and configure a Mordred calculator."""
    calc = Calculator(descriptors.all, ignore_3D=False)
    # Register additional PathCount descriptors for flexibility
    for i in range(1, 51):
        try:
            calc.register(PathCount(order=i, pi=False))
            calc.register(PathCount(order=i, pi=True))
        except: pass
    return calc

def main():
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser(description="Reproducible Descriptor Calculator.")
    parser.add_argument('--config', help='JSON configuration file.')
    parser.add_argument('--batch_size', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42, help='Seed for 3D generation.')
    args = parser.parse_args()

    if args.config:
        with open(args.config, 'r') as f: config = json.load(f)
        in_path, out_path = config['input_path'], config['output_path']
        smiles_col = config.get('smiles_col', 'smiles')
        name_col = config.get('name_col')
    else:
        sys.exit("Please provide --config file. See README for example.")

    calc = Calculator(descriptors.all, ignore_3D=False)
    # Register additional PathCount descriptors
    for i in range(1, 51):
        try:
            calc.register(PathCount(order=i, pi=False))
            calc.register(PathCount(order=i, pi=True))
        except: pass

    print(f"Reproducible mode enabled. Seed: {args.seed}")
    
    # Load and process
    df_in = pd.read_csv(in_path)
    mols, props, errors = [], [], []
    
    for i, row in tqdm(df_in.iterrows(), total=len(df_in), desc="Pre-processing"):
        mol = Chem.MolFromSmiles(str(row[smiles_col]))
        clean_mol = preprocess_molecule(mol, seed=args.seed)
        if clean_mol:
            mols.append(clean_mol)
            props.append(row.to_dict())
        else:
            errors.append({"index": i, "smiles": str(row[smiles_col])})

    if mols:
        # Mordred calculation
        m_df = calc.pandas(mols, nproc=os.cpu_count(), quiet=False)
        m_df = m_df.apply(pd.to_numeric, errors='coerce')
        
        # Custom conjugation features
        with concurrent.futures.ProcessPoolExecutor() as exec:
            c_list = list(exec.map(calc_conjugation_features, mols))
        c_df = pd.DataFrame(c_list)
        
        final_df = pd.concat([pd.DataFrame(props), m_df, c_df], axis=1)
        final_df.to_csv(out_path, index=False)
        print(f"Results saved to {out_path}")
    
    if errors:
        print(f"Warning: {len(errors)} molecules failed. Check errors.log.")
        pd.DataFrame(errors).to_csv("calculation_errors.log", index=False)

if __name__ == "__main__":
    main()
