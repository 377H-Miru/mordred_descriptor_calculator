# Extended Descriptor Calculator (Extended Mordred & Conjugation Descriptors) / 拡張記述子計算スクリプト

[English](#english) | [日本語](#japanese)

---

<a name="english"></a>
## English

### Overview
This script calculates extended molecular descriptors from SDF or CSV/TSV/SMILES files. It computes over 1,900 descriptors including all standard `mordred` descriptors (extended with PathCounts up to 50th order) and custom-developed physicochemical descriptors specialized for **$\pi$-conjugation systems**.

The tool is built for high-performance and robust processing of large-scale datasets (up to millions of molecules) using parallel processing and incremental batch saving.

### Key Features
- **Extended Mordred Descriptors**: Includes all standard descriptors plus extended PathCount (MPC/piPC) descriptors from 1st up to **50th order**.
- **Custom Conjugation Descriptors**: Specialized features for evaluating $\pi$-conjugation quality (BLA, Graph Energy, etc.).
- **Robustness for Large Datasets**: Chunked processing and incremental saving to prevent memory issues and data loss.
- **Interactive Mode**: Includes a built-in file browser for easy parameter setup.

### Installation
```bash
pip install rdkit-pypi mordred pandas numpy networkx tqdm
```

### Usage

#### 1. Interactive Mode (Recommended)
Simply run the script to start the interactive file and directory selector:
```bash
python calculate_descriptors.py
```

#### 2. Batch Mode (Using Config)
You can provide a JSON configuration file for non-interactive use:
```bash
python calculate_descriptors.py --config job_config.json
```

**Example `job_config.json`**:
```json
{
  "input_path": "your_molecules.sdf",
  "output_path": "./output_directory",
  "smiles_col": "SMILES",
  "name_col": "ID"
}
```

---

<a name="japanese"></a>
## 日本語

### 概要
本スクリプトは、SDFまたはCSV/TSV/SMILESファイルから拡張された分子記述子を計算します。標準の `mordred` 記述子（50次までのPathCount拡張を含む）と、**$\pi$共役系**に特化したカスタム物理化学記述子を含む、1,900種類以上の特徴量を算出します。

### 主な機能
- **拡張Mordred記述子**: PathCount記述子（MPC/piPC）を**1次から50次**まで拡張して計算。
- **カスタム共役系記述子**: $\pi$共役系の品質（BLA, グラフエネルギー等）を評価する独自特徴量。
- **大規模データへの対応**: 分割処理（チャンク）と逐次保存により、数百万件規模のデータに対しても安定して動作。
- **対話型モード**: ファイルブラウザ機能を内蔵しており、引数なしで簡単に実行可能。

### インストール方法
```bash
pip install rdkit-pypi mordred pandas numpy networkx tqdm
```

### 使い方

#### 1. 対話型モード (推奨)
引数なしで実行すると、入力ファイルと出力先を対話形式で選択できます：
```bash
python calculate_descriptors.py
```

#### 2. バッチモード (設定ファイルを使用)
非対話形式で実行する場合は、JSON形式の設定ファイルを指定します：
```bash
python calculate_descriptors.py --config job_config.json
```
