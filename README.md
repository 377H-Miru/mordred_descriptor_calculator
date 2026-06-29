# Mordred Descriptor Calculator CLI (`mordred-desc`)

![CI](https://github.com/377H-Miru/mordred_descriptor_calculator/actions/workflows/python-package.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

分子SMILESデータから **Mordred 記述子**（2D/3D）および**カスタム共役系記述子（Conjugation features）**を計算・抽出する公開グレードのコマンドラインツール（CLI）です。

## 動作要件 (Requirements)
- **Python 3.10 推奨 (Python 3.10 recommended)**
- (※ Python 3.11 以上は未検証・実験的サポートです)

---

## 主な機能・特徴
- **柔軟な計算モード**: デフォルトは計算が高速で安定した 2D モード。必要に応じて `--include-3d` や `--include-conjugation` を指定可能（`--only-2d` と `--include-3d` は排他的オプション）。
- **構造化エラーハンドリングと透過性**: SMILESパース、分子標準化、脱塩、3D座標生成、力場最適化、Mordred計算の各段階（`stage`）で発生した例外を記録し、`.errors.csv` / `.errors.tsv` へ出力。共役系計算中の例外も `Conjugation_Error` 列へ明示的に記録され、正常な非共役分子と区別されます。
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
  "include_conjugation": true,
  "workers": 2,
  "standardize": true,
  "desalt": true,
  "overwrite": true
}
```
実行:
```bash
mordred-desc --config job_config.json
```

#### Config JSON で指定可能なパラメータ一覧:
CLIの全オプションに対応しています。CLI引数と併用された場合は CLI引数が優先されます。
- `input_path` / `input`, `output_path` / `output`, `smiles_col`, `id_col` / `name_col`
- `seed`, `workers`, `chunksize`, `no_optimize`, `standardize`, `desalt`
- `include_3d`, `include_conjugation`, `mordred_only`, `keep_input_cols`, `minimal_output`, `overwrite`, `output_format`

---

## 計算モードと 2D/3D 記述子について
- **2D モード（デフォルト）**: 2D分子構造のみに基づく記述子を計算します。高速で失敗リスクが極めて低いモードです（明示する場合は `--only-2d`）。
- **`--include-3d`**: ETKDGv3 アルゴリズムにより 3D 座標を生成・力場最適化し、Mordred 3D 記述子を算出します（`--only-2d` と同時指定不可）。
- **`--include-conjugation`**: $\pi$ 共役系の網羅的探索、共役系原子数、共役長、BLA（Bond Length Alternation）、グラフエネルギーなどのカスタム共役系記述子を追加算出します。

### 3D構造生成の再現性
3D座標計算を行う際、`--seed 42`（デフォルト）により乱数シードが固定されるため、同一環境・同一RDKitバージョン・同一seedでは再現性の高い 3D 構造および記述子が得られます。

---

## `--workers` と計算パフォーマンスについて
- **`--workers 1`（デフォルト）**: 最も安全で推奨される動作モードです。メモリ消費を最小限に抑えます。
- **`--workers N`（N > 1）**: 前処理・Mordred計算・共役系計算の各フェーズで段階的に並列化が行われます。
  - 大規模データでは CPU およびメモリ消費量が増加するため、共有サーバーやリソース制限のある環境では小さな値から試行してください。

---

## 出力仕様およびエラーログ

### 正常出力 (`output.csv` / `output.tsv`)
- デフォルト（`--keep-input-cols` 有効）では、元の入力列すべてに加え、`canonical_smiles`（脱塩・標準化適用後の正規化SMILES）および計算された記述子列が出力されます。
- 元の入力列を出力しない場合は `--no-keep-input-cols` または `--minimal-output` を指定します（`ID`, `canonical_smiles`, および記述子列のみ出力）。

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

## サードパーティ依存関係 (Third-party dependencies)
本ツールは、SMILESのパース、分子前処理、座標計算に RDKit、分子記述子計算に Mordred を利用しています。
RDKit は BSD 系のオープンソースライセンスで配布されています。

---

## ライセンス (License)
本プロジェクトは [MIT License](LICENSE) のもとで公開されています。

---

## 引用・参照情報 (Citation)
本ツールで算出される Mordred 記述子を利用して研究成果を公表される際は、Mordred の原著論文を引用してください。

> Moriwaki, H., Tian, YS., Kawashita, N. et al. Mordred: a molecular descriptor calculator. *J Cheminform* **10**, 4 (2018). https://doi.org/10.1186/s13321-018-0258-y
