import pytest
import subprocess
import os
import sys
import json
import pandas as pd

CMD = [sys.executable, "-m", "mordred_descriptor_calculator.cli"]

def test_cli_help():
    res = subprocess.run(CMD + ["--help"], capture_output=True, text=True)
    assert res.returncode == 0
    assert "Mordred Descriptor Calculator" in res.stdout

def test_cli_execution_and_overwrite(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    
    df = pd.DataFrame({"smiles": ["CCO", "INVALID", "c1ccccc1"], "id": ["M1", "M2", "M3"]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--smiles-col", "smiles", "--id-col", "id", "--only-2d", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert os.path.exists(output_csv)
    
    err_csv = tmp_path / "output.csv.errors.csv"
    assert os.path.exists(err_csv)
    
    out_df = pd.read_csv(output_csv)
    assert len(out_df) == 2  # 2 valid SMILES
    
    err_df = pd.read_csv(err_csv)
    assert len(err_df) == 1  # 1 invalid SMILES
    assert err_df.iloc[0]["stage"] == "parse"

def test_cli_mutually_exclusive_2d_3d():
    res = subprocess.run(CMD + ["--input", "dummy.csv", "--output", "out.csv", "--only-2d", "--include-3d"], capture_output=True, text=True)
    assert res.returncode != 0
    assert "not allowed with argument" in res.stderr or "not allowed with argument" in res.stdout

def test_cli_include_3d_no_optimize(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    df = pd.DataFrame({"smiles": ["CCO"]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--include-3d", "--no-optimize", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert os.path.exists(output_csv)

def test_cli_include_conjugation(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    df = pd.DataFrame({"smiles": ["c1ccccc1"]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--include-conjugation", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    out_df = pd.read_csv(output_csv)
    assert "Conjugation_Count" in out_df.columns
    assert out_df.iloc[0]["Conjugation_Count"] >= 1

def test_cli_config_and_tsv_format(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_tsv = tmp_path / "output.tsv"
    config_json = tmp_path / "config.json"
    
    df = pd.DataFrame({"smiles": ["CCO"], "mol_id": ["A1"]})
    df.to_csv(input_csv, index=False)
    
    cfg = {
        "input_path": str(input_csv),
        "output_path": str(output_tsv),
        "smiles_col": "smiles",
        "id_col": "mol_id",
        "workers": 1,
        "standardize": True,
        "desalt": True
    }
    with open(config_json, "w", encoding="utf-8") as f:
        json.dump(cfg, f)
        
    cmd = CMD + ["--config", str(config_json), "--output-format", "tsv", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert os.path.exists(output_tsv)
    out_df = pd.read_csv(output_tsv, sep="\t")
    assert len(out_df) == 1

def test_cli_minimal_output_and_no_keep_input_cols(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    
    df = pd.DataFrame({"smiles": ["CCO", "CCC"], "extra_col": ["val1", "val2"]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--no-keep-input-cols", "--workers", "2", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    out_df = pd.read_csv(output_csv)
    assert "extra_col" not in out_df.columns
    assert "ID" in out_df.columns
    assert "canonical_smiles" in out_df.columns
