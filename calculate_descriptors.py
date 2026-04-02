import os
import sys
import warnings
import argparse
import json
import concurrent.futures
import math

# Monkey-patch for mordred compatibility with modern numpy (>= 1.25)
import numpy as np
if not hasattr(np, 'product'):
    np.product = np.prod

import networkx as nx
import pandas as pd
from mordred import Calculator, descriptors
from mordred.error import DuplicatedDescriptorName
from mordred.PathCount import PathCount
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

# =============================================================================
# 1. Interactive UI Functions
# =============================================================================
def interactive_selector(start_path, prompt, select_file=True, allowed_exts=None):
    """Interactively select a file or directory from the terminal."""
    if allowed_exts is None:
        allowed_exts = ['.sdf', '.sd', '.csv', '.tsv', '.smi', '.txt']
    current_path = os.path.abspath(start_path)

    while True:
        print("\n" + "="*50 + f"\n{prompt}\nCurrent Path: {current_path}\n" + "="*50)
        try:
            items = sorted(os.listdir(current_path))
        except OSError as e:
            print(f"Error: {e}"); current_path = os.path.dirname(current_path); continue

        options = {'0': ('.. (Parent Directory)', os.path.join(current_path, '..'))}
        if not select_file:
            options['1'] = ('./ (Select this directory)', current_path)

        dir_list, file_list = [], []
        for item in items:
            _, ext = os.path.splitext(item)
            ext = ext.lower()
            full_path = os.path.join(current_path, item)
            if os.path.isdir(full_path):
                dir_list.append((f"[DIR] {item}", full_path))
            elif select_file and ext in allowed_exts:
                file_list.append((f"[{ext.upper()[1:]}] {item}", full_path))
        
        display_list = sorted(dir_list) + sorted(file_list)
        for i, (name, _) in options.items(): print(f"[{i}] {name}")
        for i, (name, _) in enumerate(display_list, start=len(options)): print(f"[{i}] {name}")

        try:
            choice = input("Enter number to select, or 'q' to quit: ").strip()
            if choice.lower() == 'q': sys.exit("Selection cancelled.")
            choice_idx = int(choice)

            if str(choice_idx) in options:
                selected_path = options[str(choice_idx)][1]
                if str(choice_idx) == '1' and not select_file: return os.path.abspath(selected_path)
                current_path = os.path.abspath(selected_path)
            elif choice_idx >= len(options):
                selected_path = display_list[choice_idx - len(options)][1]
                if os.path.isdir(selected_path): current_path = selected_path
                else: return os.path.abspath(selected_path)
            else: print("Invalid choice.")
        except (ValueError, IndexError): print("Invalid input.")
        except Exception as e: print(f"An error occurred: {e}")

def select_column_interactive(columns, prompt_message):
    """Interactively select a column from a list of column names."""
    while True:
        print(f"\n{prompt_message}")
        for i, col in enumerate(columns): print(f"[{i}] {col}")
        try:
            choice = input("Enter the number for the column: ").strip()
            return columns[int(choice)]
        except (ValueError, IndexError): print("Invalid input. Please enter a valid number.")

# =============================================================================
# 2. Preprocessing and Standardization
# =============================================================================
def preprocess_molecule(mol):
    """Performs desalting and structure standardization."""
    if mol is None: return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        mol = Chem.MolFromSmiles(canon_smiles)
        return mol
    except Exception: return None

# =============================================================================
# 2. Descriptor Calculation Functions
# =============================================================================
def calc_conjugation_features(mol):
    """Calculates features related to conjugation systems."""
    default_features = {"Conjugation_Count": 0, "Conjugation_MaxAtomCount": 0, "Conjugation_MaxLength": 0, "Conjugation_BLA": np.nan, "Conjugation_GraphEnergy": 0.0}
    if mol is None: return default_features
    try:
        conjugated_bonds = [bond for bond in mol.GetBonds() if bond.GetIsConjugated()]
        if not conjugated_bonds: return default_features
        has_3d = mol.GetNumConformers() > 0
        conf = mol.GetConformer() if has_3d else None
        graph = nx.Graph()
        for bond in conjugated_bonds:
            start_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            weight = conf.GetAtomPosition(start_idx).Distance(conf.GetAtomPosition(end_idx)) if has_3d else 1.0
            graph.add_edge(start_idx, end_idx, weight=weight)
        conjugated_systems = list(nx.connected_components(graph))
        if not conjugated_systems: return default_features
        largest_cs = max(conjugated_systems, key=len)
        subgraph = graph.subgraph(largest_cs)
        max_length = 0
        if subgraph.number_of_nodes() > 1:
            try: max_length = nx.diameter(subgraph)
            except: pass
        bla = np.nan
        if has_3d and subgraph.number_of_edges() > 0:
            lengths = list(nx.get_edge_attributes(subgraph, 'weight').values())
            if lengths: bla = np.std(lengths)
        adj_matrix = nx.to_numpy_array(subgraph, weight=None)
        eigenvalues = np.linalg.eigvalsh(adj_matrix)
        graph_energy = np.sum(np.abs(eigenvalues))
        return {"Conjugation_Count": len(conjugated_systems), "Conjugation_MaxAtomCount": len(largest_cs), "Conjugation_MaxLength": max_length, "Conjugation_BLA": bla, "Conjugation_GraphEnergy": graph_energy}
    except Exception: return default_features

def setup_mordred_calculator():
    """Creates and configures a Mordred calculator."""
    calc = Calculator(descriptors.all, ignore_3D=False)
    for i in range(1, 51):
        try: calc.register(PathCount(order=i, pi=False)); calc.register(PathCount(order=i, pi=True))
        except DuplicatedDescriptorName: pass
    return calc

# =============================================================================
# 3. Batch Loading and Processing
# =============================================================================
def get_molecule_batches(file_path, batch_size, smiles_col=None, name_col=None, is_interactive=False):
    """Generator that yields batches of (mols, properties_df) from SDF or CSV/TSV."""
    file_ext = os.path.splitext(file_path)[1].lower()
    
    if file_ext in ['.sdf', '.sd']:
        suppl = Chem.SDMolSupplier(file_path, removeHs=False)
        current_mols, current_props = [], []
        for i, mol in enumerate(suppl):
            if mol is not None:
                clean_mol = preprocess_molecule(mol)
                if clean_mol:
                    props = {p: mol.GetProp(p) for p in mol.GetPropNames()}
                    props['Molecule_Name'] = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Mol_{i+1}"
                    if clean_mol.GetNumConformers() == 0:
                        clean_mol = Chem.AddHs(clean_mol)
                        params = AllChem.ETKDG()
                        params.maxAttempts = 100
                        AllChem.EmbedMolecule(clean_mol, params)
                    clean_mol.SetProp("_Name", props['Molecule_Name'])
                    current_mols.append(clean_mol); current_props.append(props)
            
            if len(current_mols) >= batch_size:
                yield current_mols, pd.DataFrame(current_props)
                current_mols, current_props = [], []
        if current_mols: yield current_mols, pd.DataFrame(current_props)

    elif file_ext in ['.csv', '.tsv', '.smi', '.txt']:
        delimiter = '\t' if file_ext in ['.tsv', '.smi'] else ','
        
        # In interactive mode, select columns if not provided
        if is_interactive and not smiles_col:
            sample_df = pd.read_csv(file_path, delimiter=delimiter, nrows=1)
            smiles_col = select_column_interactive(sample_df.columns, "Select the column containing SMILES strings:")
            if input("Do you want to select a column for molecule names/IDs? (y/n): ").lower() == 'y':
                name_col = select_column_interactive(sample_df.columns, "Select the column for molecule names/IDs:")

        chunks = pd.read_csv(file_path, delimiter=delimiter, chunksize=batch_size, encoding='utf-8')
        for chunk_df in chunks:
            current_mols, current_props = [], []
            if smiles_col not in chunk_df.columns:
                sys.exit(f"Error: SMILES column '{smiles_col}' not found.")
            
            for i, row in chunk_df.iterrows():
                smiles = str(row[smiles_col])
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    clean_mol = preprocess_molecule(mol)
                    if clean_mol:
                        clean_mol = Chem.AddHs(clean_mol)
                        params = AllChem.ETKDG()
                        params.maxAttempts = 100
                        AllChem.EmbedMolecule(clean_mol, params)
                        name = str(row[name_col]) if name_col and name_col in row and pd.notna(row[name_col]) else f"Mol_{i+1}"
                        clean_mol.SetProp("_Name", name)
                        current_mols.append(clean_mol)
                        props = {'Molecule_Name': name}
                        props.update({col: row[col] for col in chunk_df.columns if col != smiles_col and col != name_col})
                        current_props.append(props)
            if current_mols: yield current_mols, pd.DataFrame(current_props)

# =============================================================================
# 4. Main Execution
# =============================================================================
def main():
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser(description="Robust Batch Descriptor Calculator.")
    parser.add_argument('--config', type=str, help='Path to a JSON config file.')
    parser.add_argument('--batch_size', type=int, default=1000, help='Number of molecules per batch.')
    args = parser.parse_args()

    smiles_col, name_col = None, None
    is_interactive = False

    if args.config:
        with open(args.config, 'r', encoding='utf-8') as f: config = json.load(f)
        input_path, output_dir = os.path.abspath(config['input_path']), os.path.abspath(config['output_path'])
        smiles_col, name_col = config.get('smiles_col'), config.get('name_col')
    else:
        print("Welcome to the Extended Descriptor Calculator (Interactive Mode).")
        input_path = interactive_selector('.', "Select the input SDF, CSV, or TSV file.", select_file=True)
        output_dir = interactive_selector(os.path.dirname(input_path), "Select the output directory.", select_file=False)
        is_interactive = True

    output_csv = os.path.join(output_dir, os.path.splitext(os.path.basename(input_path))[0] + "_descriptors.csv")
    calc = setup_mordred_calculator()
    n_cores = os.cpu_count() or 1
    
    print(f"\nStarting Processing. Batch size: {args.batch_size}")
    print(f"Input: {input_path}\nOutput: {output_csv}\n")

    first_batch = True
    total_processed = 0

    # Process in batches
    batch_gen = get_molecule_batches(input_path, args.batch_size, smiles_col, name_col, is_interactive)
    
    with open(output_csv, 'w', encoding='utf-8', newline='') as f_out:
        for mols, props_df in batch_gen:
            # 1. Mordred
            m_df = calc.pandas(mols, nproc=n_cores, quiet=True)
            m_df = m_df.apply(pd.to_numeric, errors='coerce')

            # 2. Custom Conjugation (Parallel)
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
                c_list = list(executor.map(calc_conjugation_features, mols))
            c_df = pd.DataFrame(c_list)

            # 3. Combine
            props_df.reset_index(drop=True, inplace=True)
            m_df.reset_index(drop=True, inplace=True)
            c_df.reset_index(drop=True, inplace=True)
            final_df = pd.concat([props_df, m_df, c_df], axis=1)

            # Log transformation
            for col in ['Conjugation_Count', 'Conjugation_MaxAtomCount', 'Conjugation_MaxLength']:
                if col in final_df.columns:
                    final_df[f"{col}_log"] = np.log1p(final_df[col].fillna(0))

            # 4. Save (Append mode)
            final_df.to_csv(f_out, header=first_batch, index=False)
            
            total_processed += len(mols)
            first_batch = False
            print(f"Processed and saved {total_processed} molecules...")

    print("\n" + "="*50)
    print(f"Done! All results saved to: {output_csv}")
    print("="*50)

if __name__ == "__main__":
    main()
