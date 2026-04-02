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
- **Custom Conjugation Descriptors**:
    - `Conjugation_Count`: Number of independent $\pi$-conjugation systems.
    - `Conjugation_MaxAtomCount`: Number of atoms in the largest conjugation system.
    - `Conjugation_MaxLength`: Diameter (longest path) of the largest conjugation system.
    - `Conjugation_BLA`: Bond Length Alternation (std deviation of bond lengths) for the largest system (3D required).
    - `Conjugation_GraphEnergy`: Graph energy of the largest system.
- **Robustness for Large Datasets**:
    - **Chunked Processing**: Maintains low memory footprint by processing molecules in batches (default 1,000).
    - **Incremental Saving**: Automatically saves results to CSV after each batch. No data loss even if interrupted.
- **Automatic Preprocessing**: Performs desalting (keeping the largest fragment) and standardization using RDKit's `MolStandardize`.
- **Parallel Computing**: Fully utilizes all CPU cores for both Mordred and custom descriptor calculations.

### Installation
1. **Create Conda Environment**:
   ```bash
   conda create -n descriptor_env python=3.11 -y
   conda activate descriptor_env
   ```
2. **Install RDKit**:
   ```bash
   conda install -c conda-forge rdkit -y
   ```
3. **Install Dependencies**:
   ```bash
   pip install mordred pandas numpy networkx tqdm
   ```

### Usage
Run the script with a JSON configuration file:
```bash
python calculate_descriptors.py --config job_config.json --batch_size 1000
```

**Example `job_config.json`**:
```json
{
  "input_path": "your_molecules.sdf",
  "output_path": "./output_directory",
  "smiles_col": "SMILES",  // Only for CSV/TSV
  "name_col": "ID"         // Only for CSV/TSV
}
```

---

<a name="japanese"></a>
## 日本語

### 概要
本スクリプトは、SDFまたはCSV/TSV/SMILESファイルから拡張された分子記述子を計算します。標準の `mordred` 記述子（50次までのPathCount拡張を含む）と、**$\pi$共役系**に特化したカスタム物理化学記述子を含む、1,900種類以上の特徴量を算出します。

並列処理とチャンク単位の逐次保存を採用しており、数百万件規模の大規模なデータセットに対しても高いパフォーマンスと安定性を発揮します。

### 主な機能
- **拡張Mordred記述子**: 全ての標準記述子に加え、PathCount記述子（MPC/piPC）を**1次から50次**まで拡張して計算します。
- **カスタム共役系記述子**:
    - `Conjugation_Count`: 分子内の独立した共役系の数。
    - `Conjugation_MaxAtomCount`: 最大共役系に含まれる原子の数。
    - `Conjugation_MaxLength`: 最大共役系の直径（最長パス）。
    - `Conjugation_BLA`: 最大共役系の結合長交互（Bond Length Alternation）。3D構造が必要です。
    - `Conjugation_GraphEnergy`: 最大共役系のグラフエネルギー（電子的な安定性の指標）。
- **大規模データへの対応と堅牢性**:
    - **チャンク処理**: データを分割して処理（デフォルト1,000件）することで、メモリ不足を防ぎます。
    - **逐次保存 (Append Mode)**: 各バッチの計算終了ごとにCSVへ追記保存するため、途中で中断してもデータが失われません。
- **自動前処理**: RDKitの `MolStandardize` を使用し、脱塩（最大フラグメントの抽出）と構造の標準化を自動で行います。
- **並列処理**: Mordredおよびカスタム記述子の両方で、全CPUコアを活用した並列計算を行います。

### インストール方法
1. **Conda環境の作成**:
   ```bash
   conda create -n descriptor_env python=3.11 -y
   conda activate descriptor_env
   ```
2. **RDKitのインストール**:
   ```bash
   conda install -c conda-forge rdkit -y
   ```
3. **依存ライブラリのインストール**:
   ```bash
   pip install mordred pandas numpy networkx tqdm
   ```

### 使用方法
JSON形式の設定ファイルを指定して実行します。
```bash
python calculate_descriptors.py --config job_config.json --batch_size 1000
```

**`job_config.json` の例**:
```json
{
  "input_path": "your_molecules.sdf",
  "output_path": "./output_directory",
  "smiles_col": "SMILES",  // CSV/TSVの場合のみ
  "name_col": "ID"         // CSV/TSVの場合のみ
}
```

---

## Q&A / よくある質問

**Q: Computational time is too long for large molecules. / 計算に時間がかかります。**
**A:** You can disable 3D descriptors by modifying `ignore_3D=True` in `setup_mordred_calculator` to speed up the process. / 3D記述子が不要な場合は、`setup_mordred_calculator` 内で `ignore_3D=True` に設定することで高速化が可能です。

**Q: `Conjugation_BLA` values are NaN. / Conjugation_BLA が NaN になります。**
**A:** This descriptor requires 3D coordinates. Ensure your SDF has 3D info or wait for the script's automatic 3D generation from SMILES. / この記述子には3D座標が必要です。SDFに3D情報が含まれているか確認してください。SMILES入力の場合は自動生成を試みますが、生成に失敗した場合はNaNとなります。
