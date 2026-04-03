# mordred_descriptor_calculator

Mordred記述子に加え、$\pi$共役系に特化したカスタム記述子を算出する堅牢なツール。

## 特徴
- **再現性**: 3D構造生成のシード値を固定（デフォルト: 42）。
- **堅牢性**: 不正な構造や3D Embed失敗を安全にスキップ。
- **専門性**: $\pi$共役系の品質（Conjugation Features）を独自に算出。

## 使い方 (Usage)

### 1. コマンドライン実行 (CLI Mode)
JSON設定ファイルを作成し、引数として渡します。
```bash
python calculate_descriptors.py --config job_config.json --seed 42
```

#### job_config.json の例:
```json
{
  "input_path": "compounds.csv",
  "output_path": "results_descriptors.csv",
  "smiles_col": "smiles",
  "name_col": "compound_id"
}
```

### 2. インタラクティブモード
引数なしで実行すると、フォルダをブラウズしてファイルを選択できます。
```bash
python calculate_descriptors.py
```

## インストール
```bash
pip install -r requirements.txt
```

## 算出されるカスタム記述子
- `Conjugation_Count`: 共役系の数
- `Conjugation_MaxAtomCount`: 最大の共役系のサイズ
- `Conjugation_MaxLength`: 共役系の最長経路
- `Conjugation_BLA`: 結合長交互作用 (Bond Length Alternation)
- `Conjugation_GraphEnergy`: グラフエネルギー
