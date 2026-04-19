import os
import sys
import warnings
import argparse
import json
import concurrent.futures
import numpy as np

# --- ULTIMATE MONKEY-PATCH FOR MORDRED & MODERN NUMPY ---
if not hasattr(np, 'product'): np.product = np.prod
import numpy.core.fromnumeric
sys.modules['numpy.product'] = np.prod
# --------------------------------------------------------

import networkx as nx
import pandas as pd
from mordred import Calculator, descriptors
from mordred.PathCount import PathCount
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

def preprocess_worker(args):
    """Worker function for parallel preprocessing with deterministic seed."""
    smiles, props, seed, optimize = args
    if pd.isna(smiles) or str(smiles).strip() == "" or str(smiles) == "nan":
        return None, None, {"smiles": "Empty", "reason": "Data Error: SMILES is empty"}
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None: return None, None, {"smiles": smiles, "reason": "RDKit Parse Error"}
        mol = rdMolStandardize.Cleanup(mol)
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(mol, params) == -1: return None, None, {"smiles": smiles, "reason": "3D Embedding Failed"}
        if optimize:
            try:
                if AllChem.MMFFOptimizeMolecule(mol, maxIters=200) == -1:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except: pass
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
    parser = argparse.ArgumentParser(description="Professional & Reproducible Descriptor Calculator.")
    parser.add_argument('--config', help='JSON configuration file.')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--no-optimize', action='store_true')
    parser.add_argument('--chunksize', type=int, default=1000)
    args = parser.parse_args()

    if not args.config:
        print("\n--- Mordred Calculator (Interactive Mode) ---")
        in_path = input("Enter input CSV/TSV path > ").strip()
        if not os.path.exists(in_path): sys.exit("File not found.")
        out_path = input("Enter output CSV path [default: results.csv] > ").strip() or "results.csv"
        smiles_col = input("Enter SMILES column name [default: smiles] > ").strip() or "smiles"
    else:
        with open(args.config, 'r', encoding='utf-8') as f: config = json.load(f)
        in_path, out_path = config['input_path'], config['output_path']
        smiles_col = config.get('smiles_col', 'smiles')

    calc = setup_mordred_calculator()
    global_count = 0
    first_chunk = True
    n_workers = os.cpu_count() or 1

    print(f"Processing in chunks of {args.chunksize} using {n_workers} workers...")
    
    # Process in chunks to handle huge datasets
    for chunk in pd.read_csv(in_path, sep=None, engine='python', chunksize=args.chunksize):
        tasks = []
        for i, row in enumerate(chunk.to_dict('records')):
            tasks.append((row[smiles_col], row, args.seed + global_count + i, not args.no_optimize))
        
        mols, props, errors = [], [], []
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(tqdm(executor.map(preprocess_worker, tasks), total=len(tasks), desc=f"Chunk {global_count//args.chunksize + 1} Prep", leave=False))

        for mol, p, err in results:
            if mol:
                mols.append(mol); props.append(p)
            else:
                errors.append(err)

        if mols:
            m_df = calc.pandas(mols, nproc=n_workers, quiet=True)
            m_df = m_df.apply(pd.to_numeric, errors='coerce')
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
                c_list = list(executor.map(calc_conjugation_features, mols))
            c_df = pd.DataFrame(c_list)
            
            final_df = pd.concat([pd.DataFrame(props).reset_index(drop=True), m_df.reset_index(drop=True), c_df.reset_index(drop=True)], axis=1)
            final_df.to_csv(out_path, index=False, header=first_chunk, mode='a' if not first_chunk else 'w', encoding='utf-8-sig')
            first_chunk = False
        
        if errors:
            err_file = out_path + ".errors.log"
            pd.DataFrame(errors).to_csv(err_file, index=False, mode='a', header=not os.path.exists(err_file))
            
        global_count += len(chunk)

    print(f"\nDone. Results saved to {out_path}")

if __name__ == "__main__":
    main()
