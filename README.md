# Mordred Descriptor Calculator CLI (`mordred-desc`)

![CI](https://github.com/377H-Miru/mordred_descriptor_calculator/actions/workflows/python-package.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

分子SMILESデータから **Mordred 記述子**（2D/3D）および**カスタム共役系記述子（Conjugation features）**を高速かつ堅牢に計算・算出する公開グレードのコマンドラインツール（CLI）です。

## 主な機能・特徴
- **柔軟な計算モード**: デフォルトは計算が高速で安定した `--only-2d` モード。必要に応じて `--include-3d` や `--include-conjugation` を指定可能。
- **堅牢な多段階エラーハンドリング**: SMILESパース、分子標準化、脱塩、3D座標生成（ETKDGv3）、構造最適化（MMFF/UFF）の各段階で発生したエラーを詳細に分類し、`.errors.csv` / `.errors.tsv` へ出力。
- **安全な並列処理管理**: 前処理フェーズと Mordred 計算フェーズでプロセス競合（Nested Multiprocessing）が発生しないよう段階的にワーカーを制御。
- **モダンなPython互換性**: 近年の NumPy バージョンで非推奨・削除された `np.product` 等のエイリアスに対する自動互換性シムを内蔵。

---

## インストール (Installation)

```bash
git clone https://github.com/377H-Miru/mordred_descriptor_calculator.git
cd mordred_descriptor_calculator
pip install -e .
```

ヘルプメッセージの表示:
```bash
mordred-desc --help
```

---

## Quick Start & 実行例

### 1. 最小入力CSV例 (`sample.csv`)
```csv
smiles,compound_id
CCO,MOL001
c1ccccc1,MOL002
CO.[Na+].[Cl-],MOL003
INVALID_SMILES,MOL004
```

### 2. コマンドライン実行例 (CLI Mode)
```bash
mordred-desc \
  --input sample.csv \
  --output result.csv \
  --smiles-col smiles \
  --id-col compound_id \
  --only-2d \
  --overwrite
```

### 3. 設定ファイル実行例 (Config JSON Mode)
`job_config.json`:
```json
{
  "input_path": "sample.csv",
  "output_path": "result.csv",
  "smiles_col": "smiles",
  "id_col": "compound_id",
  "include_3d": false,
  "include_conjugation": true
}
```
実行:
```bash
mordred-desc --config job_config.json --overwrite
```

---

## 計算モードと 2D/3D 記述子について
- **`--only-2d`（デフォルト推薦）**: 2D分子構造のみに基づく記述子を計算します。高速で失敗リスクが極めて低いモードです。
- **`--include-3d`**: ETKDGv3 アルゴリズムにより 3D 座標を生成・力場最適化し、Mordred 3D 記述子を算出します。
- **`--include-conjugation`**: $\pi$ 共役系の網羅的探索、共役系原子数、共役長、BLA（Bond Length Alternation）、グラフエネルギーなどのカスタム共役系記述子を追加算出します。

### 3D構造生成の再現性
3D座標計算を行う際、`--seed 42`（デフォルト）により乱数シードが固定されるため、同一入力に対して完全に再現可能な 3D 構造および記述子が得られます。

---

## `--workers` と計算パフォーマンス
- `--workers 1`（デフォルト）: シングルプロセスで逐次処理を行います。メモリ消費が最も少なく安全です。
- `--workers N`（N > 1）: `ProcessPoolExecutor` により前処理および Mordred 計算をマルチプロセス化します。前処理完了後に Mordred の並列計算を開始する設計となっており、CPUsの過剰消費やデッドロックを防ぎます。

---

## 出力仕様およびエラーログ

### 正常出力 (`output.csv` / `output.tsv`)
- デフォルト（`--keep-input-cols` 有効）では、元の入力列すべてに加え、`canonical_smiles`（脱塩・標準化適用後の正規化SMILES）および計算された記述子列が出力されます。
- `--minimal-output` を指定すると、`ID`, `canonical_smiles`, および記述子列のみに削減されます。

### エラーログ (`output.csv.errors.csv` / `output.tsv.errors.tsv`)
計算過程でエラーが発生した分子は正常出力から除外され、以下の構造化フォーマットで記録されます。
- `row_index`: 入力データの行番号
- `ID`: 化合物ID/名称
- `input_smiles`: 元の入力SMILES
- `stage`: 失敗段階（`parse`, `standardize`, `desalt`, `add_hydrogen`, `embed_3d`, `optimize_3d`, `mordred`, `conjugation`）
- `error_type`: エラー種別
- `error_message`: 詳細エラーメッセージ

---

## 安全な上書き挙動 (`--overwrite`)
出力ファイルおよびエラーログファイルが既に存在する場合、意図しないデータ消失を防ぐためデフォルトではエラー終了します。既存ファイルを上書きする場合は、必ず `--overwrite` オプションを付与してください。

---

## ライセンス (License)
本プロジェクトは [MIT License](LICENSE) のもとで公開されています。

---

## 引用・参照情報 (Citation)
本ツールで算出される Mordred 記述子を利用して研究成果を公表される際は、Mordred の原著論文を引用してください。

> Moriwaki, H., Tian, YS., Kawashita, N. et al. Mordred: a molecular descriptor calculator. *J Cheminform* **10**, 4 (2018). https://doi.org/10.1186/s13321-018-0258-y
