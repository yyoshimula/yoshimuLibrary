# Solar Radiation Pressure (SRP) Calculation Functions

衛星の太陽輻射圧（SRP: Solar Radiation Pressure）による力とトルクを計算するための関数群です．
`object` フォルダで読み込んだ衛星モデル（`sat`構造体）を用いて，様々なBRDF（双方向反射率分布関数）モデルに基づいた物理計算を行います．

## Functions

### Core SRP Models

- **`srpSimple.m`**
  Lambertian拡散反射と完全鏡面反射（Perfect Specular）を組み合わせた簡易モデルによるSRP計算関数です．
  - 解析的な計算が可能で高速です．
  - 入力: 衛星モデル `sat`，太陽方向ベクトル `sunB`，距離 `d`，定数構造体 `const`
  - 出力: SRPによる力 `sat.force`，トルク `sat.torque`

- **`srpAS.m`**
  Ashikhmin-Shirleyモデルを用いたSRP計算関数です．異方性反射を考慮できます．
  - 拡散反射成分は解析解を使用し，鏡面反射成分は重要サンプリング（Importance Sampling）を用いたモンテカルロ積分により数値的に求めます．
  - `nMC` 引数でサンプル数を指定可能です．

- **`srpCT.m`**
  Cook-Torranceモデルを用いたSRP計算関数です．微小面の粗さを考慮した計算が行えます．
  - 重要サンプリングを用いたモンテカルロ積分を使用します．
  - 法線分布関数（NDF）として 'Beckmann'（デフォルト）または 'Gauss' を選択可能です．

- **`srpCannon.m`**
  キャノンボールモデル（Cannonball model）を用いた簡易的なSRP加速度計算関数です．
  - 形状詳細を考慮せず，断面積対質量比（Area-to-Mass ratio）と反射係数（Cr）に基づき加速度を算出します．

### Variants (Sampling Methods)

- **`srpCTuni.m`**, **`srpASuni.m`**
  それぞれCook-Torrance, Ashikhmin-Shirleyモデルにおいて，一様分布（Uniform distribution）によるサンプリングを用いた旧実装または比較用関数です．
  - 重要サンプリングに比べて収束が遅いため，通常は `srpCT.m`, `srpAS.m` の使用が推奨されます．

### Examples & Verification

- **`exampleSRP.m`**
  基本的なSRP計算のデモンストレーションスクリプトです．
- **`verifyImpSampling.m`**
  各モデル（CT, AS）におけるモンテカルロ積分の収束性や，重要サンプリング（Importance Sampling）と一様分布サンプリング（Uniform Sampling）の比較検証を行うためのスクリプトです．

---

# Solar Radiation Pressure (SRP) Calculation Functions

A collection of functions for calculating forces and torques due to Solar Radiation Pressure (SRP) on a satellite.
Using a satellite model (`sat` structure) loaded in the `object` folder, it performs physical calculations based on various BRDF (Bidirectional Reflectance Distribution Function) models.

## Functions

### Core SRP Models

- **`srpSimple.m`**
  SRP calculation function using a simplified model combining Lambertian diffuse reflection and Perfect Specular reflection.
  - Analytical calculation is possible and fast.
  - Inputs: Satellite model `sat`, Sun direction vector `sunB`, distance `d`, constant structure `const`
  - Outputs: Force due to SRP `sat.force`, torque `sat.torque`

- **`srpAS.m`**
  SRP calculation function using the Ashikhmin-Shirley model. It can account for anisotropic reflection.
  - The diffuse reflection component uses an analytical solution, while the specular reflection component is determined numerically by Monte Carlo integration using Importance Sampling.
  - The number of samples can be specified with the `nMC` argument.

- **`srpCT.m`**
  SRP calculation function using the Cook-Torrance model. It allows for calculation considering the roughness of microfacets.
  - Uses Monte Carlo integration with Importance Sampling.
  - 'Beckmann' (default) or 'Gauss' can be selected as the Normal Distribution Function (NDF).

- **`srpCannon.m`**
  Simplified SRP acceleration calculation function using the Cannonball model.
  - Calculates acceleration based on Area-to-Mass ratio and reflection coefficient (Cr) without considering detailed shape.

### Variants (Sampling Methods)

- **`srpCTuni.m`**, **`srpASuni.m`**
  Old implementations or comparison functions for Cook-Torrance and Ashikhmin-Shirley models, respectively, using sampling from a Uniform distribution.
  - Convergence is slower compared to Importance Sampling, so `srpCT.m` and `srpAS.m` are generally recommended.

### Examples & Verification

- **`exampleSRP.m`**
  A demonstration script for basic SRP calculation.
- **`verifyImpSampling.m`**
  A script for verifying the convergence of Monte Carlo integration for each model (CT, AS) and comparing Importance Sampling with Uniform Sampling.
