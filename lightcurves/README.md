# Light Curve Calculation Functions

衛星のライトカーブ（光度曲線）を計算するための一般関数群です．
`object` フォルダで読み込んだ衛星モデル（`sat`構造体）を用いて，様々なBRDF（双方向反射率分布関数）モデルに基づいた光度計算を行います．

## Functions

### Main Calculation

- **`lc.m`**
  ライトカーブ計算のメイン関数です．指定されたBRDFモデル（Simple, Ashikhmin-Shirley, Cook-Torrance）に応じて，適切なサブ関数を呼び出し，衛星の姿勢や位置関係に基づいて見かけの等級（magnitude）を計算します．
  - 入力: 衛星モデル `sat`，姿勢クォータニオン `q`，位置ベクトル（衛星，観測者，太陽），日照フラグ `nu`
  - 出力: 等級 `m`，反射光量 `fObs`

### BRDF Models

各BRDFモデルに基づく反射光量の計算関数です．`lc.m` から呼び出されますが，単独で使用することも可能です．

- **`lcSimple.m`**
  Lambertian拡散反射と完全鏡面反射（Perfect Specular）を組み合わせた簡易モデルです．
  
- **`lcAS.m`**
  Ashikhmin-Shirleyモデルを用いた計算関数です．異方性反射（Anisotropic reflection）を考慮することができます．
  - 拡散反射成分と鏡面反射成分をそれぞれ計算します．

- **`lcCT.m`**
  Cook-Torranceモデルを用いた計算関数です．微小面の粗さ（Roughness）を考慮した物理ベースの計算が行えます．
  - 法線分布関数（NDF）として 'Beckmann'（デフォルト）または 'Gauss' を選択可能です．

### Utilities

- **`mag.m`**
  反射光量（Flux）と距離から，見かけの等級（Apparent Magnitude）を計算します．
  - 基準として太陽の等級 (-26.7) を使用しています．

- **`magInv.m`**
  見かけの等級と距離から，反射光量（Flux）を逆算します．

### Visualization

- **`visualizeBRDF.m`**
  設定したBRDFモデルの反射特性を可視化するためのスクリプトです．
  - 半球面上での反射強度の分布を3次元プロットで確認できます．
