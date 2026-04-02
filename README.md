# Extended Descriptor Calculator (Extended Mordred & Conjugation Descriptors)

[English](#english) | [日本語](#japanese)

---

<a name="english"></a>
## English

### Overview
This script calculates extended molecular descriptors from SDF or CSV/TSV/SMILES files. It computes over 1,900 descriptors, including all standard Mordred descriptors and custom features for $\pi$-conjugation systems.

### Installation
It is recommended to use a virtual environment (venv or conda).

1. **Clone or download** this repository.
2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

### Usage
#### 1. Interactive Mode (Recommended)
Run the script without arguments to start the interactive file browser:
```bash
python calculate_descriptors.py
```

#### 2. Batch Mode
```bash
python calculate_descriptors.py --config job_config.json
```

---

<a name="japanese"></a>
## 日本語

### 概要
SDFまたはCSV/TSV/SMILESファイルから、標準Mordred記述子および $\pi$共役系特化のカスタム記述子（合計1,900種類以上）を計算するツールです。

### インストール方法
仮想環境（venvまたはconda）の使用を推奨します。

1. **リポジトリをクローンまたはダウンロード**します。
2. **依存ライブラリのインストール**:
   ```bash
   pip install -r requirements.txt
   ```

### 使い方
#### 1. 対話型モード (推奨)
引数なしで実行すると、対話形式でファイルを選択できます：
```bash
python calculate_descriptors.py
```
