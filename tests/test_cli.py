import pytest
import subprocess
import os
import sys
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

    # Test missing --overwrite fails when output exists
    cmd_no_ow = CMD + ["--input", str(input_csv), "--output", str(output_csv)]
    res_no_ow = subprocess.run(cmd_no_ow, capture_output=True, text=True)
    assert res_no_ow.returncode != 0
    assert "already exists" in res_no_ow.stderr or "already exists" in res_no_ow.stdout

def test_cli_missing_smiles_col(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    
    df = pd.DataFrame({"structure": ["CCO"], "id": [1]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--smiles-col", "smiles"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode != 0
    combined = res.stdout + res.stderr
    assert "not found in input file" in combined
    assert "Available columns" in combined

def test_cli_workers_1_and_no_optimize(tmp_path):
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    
    df = pd.DataFrame({"smiles": ["CCO"]})
    df.to_csv(input_csv, index=False)
    
    cmd = CMD + ["--input", str(input_csv), "--output", str(output_csv), "--workers", "1", "--no-optimize", "--only-2d", "--overwrite"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert os.path.exists(output_csv)
