import argparse
import os
import sys
import json
import logging
import warnings
import concurrent.futures
import pandas as pd
import numpy as np

from .compat import patch_numpy_for_mordred_compat
from .core import preprocess_worker, setup_mordred_calculator
from .conjugation import calc_conjugation_features

def main():
    patch_numpy_for_mordred_compat()
    warnings.filterwarnings('ignore')
    
    parser = argparse.ArgumentParser(description="Mordred Descriptor Calculator - Professional Molecular Feature Extraction CLI")
    
    parser.add_argument('--input', help='Input CSV/TSV file path')
    parser.add_argument('--output', help='Output file path')
    parser.add_argument('--smiles-col', help='Column name containing SMILES (default: smiles)')
    parser.add_argument('--id-col', '--name-col', dest='id_col', help='Column name for compound ID/name')
    parser.add_argument('--config', help='JSON configuration file path')
    
    parser.add_argument('--seed', type=int, default=42, help='Random seed for 3D coordinate embedding (default: 42)')
    parser.add_argument('--workers', type=int, default=1, help='Number of parallel worker processes (default: 1)')
    parser.add_argument('--chunksize', type=int, default=1000, help='Chunk size for processing rows (default: 1000)')
    parser.add_argument('--no-optimize', action='store_true', help='Disable 3D geometry force field optimization')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output file and error log if they exist')
    
    # Standardize & Desalt groups
    std_group = parser.add_mutually_exclusive_group()
    std_group.add_argument('--standardize', action='store_true', dest='standardize', default=True, help='Apply RDKit MolStandardize cleanup (default: True)')
    std_group.add_argument('--no-standardize', action='store_false', dest='standardize', help='Disable standardization')
    
    desalt_group = parser.add_mutually_exclusive_group()
    desalt_group.add_argument('--desalt', action='store_true', dest='desalt', default=True, help='Remove salts / keep largest fragment (default: True)')
    desalt_group.add_argument('--no-desalt', action='store_false', dest='desalt', help='Disable desalting')
    
    parser.add_argument('--output-format', choices=['csv', 'tsv'], help='Output format (default: auto-detected from extension)')
    
    # Computation modes (Mutually exclusive 2D / 3D flags)
    mode_group = parser.add_argument_group('Computation Modes')
    dim_group = mode_group.add_mutually_exclusive_group()
    dim_group.add_argument('--only-2d', action='store_true', help='Explicitly specify 2D descriptors only without 3D embedding')
    dim_group.add_argument('--include-3d', action='store_true', help='Include 3D descriptors and 3D embedding')
    
    mode_group.add_argument('--include-conjugation', action='store_true', help='Include custom conjugation descriptors')
    mode_group.add_argument('--mordred-only', action='store_true', help='Calculate Mordred descriptors only (exclude conjugation)')
    
    # Input/Output filtering
    io_group = parser.add_argument_group('Output Filtering & Formatting')
    io_group.add_argument('--keep-input-cols', action='store_true', default=True, help='Keep all input columns in output (default: True)')
    io_group.add_argument('--minimal-output', action='store_true', help='Keep only ID, canonical_smiles, and descriptors in output')
    io_group.add_argument('--drop-all-na', action='store_true', help='Drop descriptor columns where all values are NaN')
    io_group.add_argument('--drop-constant', action='store_true', help='Drop descriptor columns with constant values')
    io_group.add_argument('--missing-value', choices=['nan', 'blank'], default='nan', help='Format missing descriptor values as NaN or blank string')
    
    # Logging
    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument('--quiet', action='store_true', help='Suppress progress bar and informative output')
    log_group.add_argument('--verbose', action='store_true', help='Enable verbose debug log output')
    
    args = parser.parse_args()
    
    # Setup logging level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    elif args.quiet:
        logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        
    # --- Config merging ---
    config = {}
    if args.config:
        if not os.path.exists(args.config):
            sys.exit(f"Error: Configuration file '{args.config}' not found.")
        try:
            with open(args.config, 'r', encoding='utf-8') as f:
                config = json.load(f)
        except Exception as e:
            sys.exit(f"Error parsing configuration JSON: {e}")
            
    input_file = args.input or config.get('input_path') or config.get('input')
    output_file = args.output or config.get('output_path') or config.get('output')
    smiles_col = args.smiles_col or config.get('smiles_col') or 'smiles'
    id_col = args.id_col or config.get('id_col') or config.get('name_col')
    
    if not input_file or not output_file:
        parser.print_help()
        sys.exit("\nError: Both --input and --output (or valid --config) are required.")

    # Validation
    if not os.path.exists(input_file):
        sys.exit(f"Error: Input file '{input_file}' does not exist.")
        
    if args.workers <= 0:
        sys.exit("Error: --workers must be > 0.")
    if args.chunksize <= 0:
        sys.exit("Error: --chunksize must be > 0.")

    # Computation mode flags resolution
    include_3d = args.include_3d or config.get('include_3d', False)
    if args.only_2d:
        include_3d = False
        
    ignore_3d = not include_3d
    include_conjugation = (args.include_conjugation or config.get('include_conjugation', False)) and not args.mordred_only
    optimize = not args.no_optimize and config.get('optimize', True)
    standardize = args.standardize
    desalt = args.desalt
    seed = args.seed
    chunksize = args.chunksize
    workers = args.workers
    overwrite = args.overwrite
    
    output_format = args.output_format
    if not output_format:
        ext = os.path.splitext(output_file)[1].lower()
        output_format = 'tsv' if ext in ['.tsv', '.txt'] else 'csv'
        
    err_file = output_file + (".errors.tsv" if output_format == 'tsv' else ".errors.csv")
    
    if os.path.exists(output_file) and not overwrite:
        sys.exit(f"Error: Output file '{output_file}' already exists. Use --overwrite to overwrite.")
    if os.path.exists(err_file) and not overwrite:
        sys.exit(f"Error: Error log file '{err_file}' already exists. Use --overwrite to overwrite.")

    if overwrite:
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(err_file):
            os.remove(err_file)

    sep = "\t" if input_file.lower().endswith(('.tsv', '.txt')) else ","
    
    # Peek columns
    try:
        peek_df = pd.read_csv(input_file, sep=sep, nrows=1, engine='python')
        available_cols = list(peek_df.columns)
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")
        
    if smiles_col not in available_cols:
        sys.exit(f"Error: SMILES column '{smiles_col}' not found in input file. Available columns: {available_cols}")
    if id_col and id_col not in available_cols:
        sys.exit(f"Error: ID column '{id_col}' not found in input file. Available columns: {available_cols}")

    calc = setup_mordred_calculator(ignore_3d=ignore_3d)
    
    global_count = 0
    first_chunk = True
    
    try:
        total_rows = sum(1 for _ in open(input_file, 'rb')) - 1
    except Exception:
        total_rows = None

    if not args.quiet:
        print(f"\nProcessing {input_file} -> {output_file}")
        print(f"Settings: 3D={include_3d}, Conjugation={include_conjugation}, Standardize={standardize}, Desalt={desalt}, Workers={workers}")

    for chunk in pd.read_csv(input_file, sep=sep, engine='python', chunksize=chunksize):
        tasks = []
        for i, row in enumerate(chunk.to_dict('records')):
            row_idx = global_count + i
            mol_id = str(row.get(id_col)) if id_col and id_col in row else f"ID_{row_idx}"
            tasks.append((row[smiles_col], row, seed + row_idx, optimize, standardize, desalt, include_3d, row_idx, mol_id))
            
        mols, canonical_smiles_list, props, errors = [], [], [], []
        
        if workers <= 1:
            prep_results = list(map(preprocess_worker, tasks))
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                prep_results = list(executor.map(preprocess_worker, tasks))
                
        for mol, can_smiles, p, err in prep_results:
            if mol is not None:
                mols.append(mol)
                canonical_smiles_list.append(can_smiles)
                props.append(p)
            else:
                errors.append(err)
                
        if mols:
            # Compute Mordred descriptors with error handling
            try:
                m_df = calc.pandas(mols, nproc=workers, quiet=args.quiet)
                m_df = m_df.apply(pd.to_numeric, errors='coerce')
            except Exception as me:
                # Log mordred stage failure for all mols in chunk
                for idx, p in enumerate(props):
                    mol_id = str(p.get(id_col)) if id_col and id_col in p else f"ID_{global_count + idx}"
                    errors.append({
                        'row_index': global_count + idx,
                        'ID': mol_id,
                        'input_smiles': str(p.get(smiles_col, '')),
                        'stage': 'mordred',
                        'error_type': type(me).__name__,
                        'error_message': str(me)
                    })
                if errors:
                    err_df = pd.DataFrame(errors)
                    err_sep = "," if output_format == 'csv' else "\t"
                    err_df.to_csv(err_file, sep=err_sep, index=False, mode='a', header=not os.path.exists(err_file), encoding='utf-8-sig')
                sys.exit(f"Error during Mordred calculation stage: {me}")

            # Compute conjugation features if enabled
            if include_conjugation:
                if workers <= 1:
                    c_list = list(map(calc_conjugation_features, mols))
                else:
                    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                        c_list = list(executor.map(calc_conjugation_features, mols))
                c_df = pd.DataFrame(c_list)
            else:
                c_df = pd.DataFrame()
                
            # Prepare properties / metadata DataFrame
            if args.minimal_output:
                meta_records = []
                for i, p in enumerate(props):
                    mol_id = str(p.get(id_col)) if id_col and id_col in p else f"ID_{global_count + i}"
                    meta_records.append({'ID': mol_id, 'canonical_smiles': canonical_smiles_list[i]})
                meta_df = pd.DataFrame(meta_records)
            else:
                meta_df = pd.DataFrame(props).reset_index(drop=True)
                meta_df['canonical_smiles'] = canonical_smiles_list
                
            dfs_to_concat = [meta_df.reset_index(drop=True), m_df.reset_index(drop=True)]
            if not c_df.empty:
                dfs_to_concat.append(c_df.reset_index(drop=True))
                
            final_df = pd.concat(dfs_to_concat, axis=1)
            
            # Post-processing filters
            if args.drop_all_na:
                final_df = final_df.dropna(how='all', axis=1)
            if args.drop_constant:
                numeric_cols = final_df.select_dtypes(include=[np.number]).columns
                constant_cols = [c for c in numeric_cols if final_df[c].nunique(dropna=False) <= 1]
                final_df = final_df.drop(columns=constant_cols)
            if args.missing_value == 'blank':
                final_df = final_df.fillna('')
                
            out_sep = "," if output_format == 'csv' else "\t"
            final_df.to_csv(output_file, sep=out_sep, index=False, header=first_chunk, mode='a' if not first_chunk else 'w', encoding='utf-8-sig')
            first_chunk = False
            
        if errors:
            err_df = pd.DataFrame(errors)
            err_sep = "," if output_format == 'csv' else "\t"
            err_df.to_csv(err_file, sep=err_sep, index=False, mode='a', header=not os.path.exists(err_file), encoding='utf-8-sig')
            
        global_count += len(chunk)

    if not args.quiet:
        print(f"\nSuccess. Results saved to {output_file}")

if __name__ == '__main__':
    main()
