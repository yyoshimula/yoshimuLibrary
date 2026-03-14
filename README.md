# yoshimuLibrary

本リポジトリは、宇宙工学および姿勢・軌道力学研究の文脈で開発された、様々なシミュレーション、推定アルゴリズム、制御モデルのためのMATLABコードとデータを含んでいます。

主なアプリケーションは以下の通りです：
- **太陽輻射圧 (SRP) モデリング** (`srp`): 物理ベースのBRDFモデル（Simple, Ashikhmin-Shirley, Cook-Torrance）を用いたSRP（力・トルク）の高精度計算。効率化のための簡易実装や補正係数テーブルも含みます。
- **ライトカーブシミュレーション** (`lightcurves`): 詳細な表面特性 (BRDF) や等級変換を組み込んだ、宇宙物体の合成ライトカーブ生成。
- **宇宙機オブジェクト処理** (`object`): Wavefront OBJファイルの読み込み・処理 (`readSC`)、宇宙機モデルの可視化、法線・面積計算、自己影 (`selfShadow`) の考慮などのツール群。
- **姿勢・軌道**: 推定・フィルタリング（Kalman Filter, UKF, CKFなど）、姿勢表現とその変換。

## 動作要件

- MATLAB R2024a 以降
- Git LFS がインストールされていること

## セットアップ

### リポジトリのクローン

```bash
git clone git@github.com:yyoshimula/yoshimuLibrary.git
cd yoshimuLibrary
git lfs install
git lfs pull
```

## ディレクトリ構成

- **000refs**: 参考文献・資料。
- **50attitudeKnock**: 姿勢制御・推定の演習スクリプト集。
- **50orbitKnock**: 軌道力学の演習スクリプト集。
- **SPICE**: NASA SPICE Toolkit 関連スクリプト。
- **attitude**: 姿勢力学・運動学、制御に関する関数群（DCM, クォータニオン, オイラー角など）。
- **conversion**: 単位変換や時間変換のユーティリティ（au to km, hms to degなど）。
- **dualQuaternions**: デュアルクォータニオンの演算・運動学。
- **environment**: 宇宙環境モデル（IGRF, 大気密度, 外乱など）。
- **examples**: 様々なシナリオの例題スクリプト（UKF, 姿勢・軌道など）。
- **geometricIntegration**: 幾何学的積分法（Lie群積分など）。
- **gpr**: ガウス過程回帰 (GPR) 関連関数。
- **hifiSRP**: 高精度な太陽輻射圧 (SRP) モデリング。
- **lightcurves**: ライトカーブの計算・解析。
- **math**: 一般的な数学関数（ルジャンドル多項式, 歪対称行列など）。
- **object**: 3Dオブジェクトの読み込み・処理（.obj/.mtlファイル, Self-shadowing）。
- **orbit**: 軌道力学関連関数（ケプラー方程式, 座標変換, EGM2008, EOPなど）。
- **orbitDetermination**: 軌道決定アルゴリズム（Gibbs, Gauss, Double Rなど）。
- **relativeOrbit**: 相対軌道力学・運動学（HCW, ROEなど）。
- **sphericalGaussian**: 球面ガウス関数 (Spherical Gaussian)。
- **srp**: 太陽輻射圧 (SRP) モデル（Analytical, Cannonballなど）。
- **sunMoon**: 太陽・月の位置計算、エフェメリス。
- **time**: 時間系変換（JD, MJD, UTC, TTなど）。
- **ukfCkf**: Unscented / Cubature Kalman Filter 実装。
- **utility**: 汎用ユーティリティ（プロット, 図保存など）。


---

# yoshimuLibrary

This repository contains MATLAB code and data for various simulations, estimation algorithms, and control models developed in the context of spacecraft engineering and orbital mechanics research.

Main applications include:
- **Solar Radiation Pressure (SRP) Modeling** (`srp`): High-fidelity calculation of SRP forces and torques using physically based BRDF models (Simple, Ashikhmin-Shirley, Cook-Torrance). Includes simplified implementations and correction coefficient tables for efficiency.
- **Light Curve Simulation** (`lightcurves`): Synthetic light curve generation for space objects incorporating detailed surface properties (BRDFs) and magnitude conversion.
- **Spacecraft Object Handling** (`object`): Tools for reading and processing Wavefront OBJ files (`readSC`), visualizing spacecraft models, and calculating geometric properties like normals, areas, and self-shadowing effects (`selfShadow`).
- **Attitude & Orbit**: Estimations and filtering (e.g., Kalman Filter, UKF, CKF), attitude representations, and their transformations.

---

## Requirements

- MATLAB R2024a or later
- Git LFS installed

## Setup

### Clone the repository

```bash
git clone git@github.com:yyoshimula/yoshimuLibrary.git
cd yoshimuLibrary
git lfs install
git lfs pull
```

## Directory Structure

- **000refs**: References and documents.
- **50attitudeKnock**: Collection of practice scripts for attitude control/estimation.
- **50orbitKnock**: Collection of practice scripts for orbital mechanics.
- **SPICE**: Scripts related to NASA SPICE Toolkit.
- **attitude**: Attitude dynamics, kinematics, and control functions (DCM, quaternions, Euler angles, etc.).
- **conversion**: Unit and time conversion utilities (au to km, hms to deg, etc.).
- **dualQuaternions**: Dual quaternion operations and kinematics.
- **environment**: Space environment models (IGRF, atmospheric density, disturbances).
- **examples**: Example scripts for various scenarios (UKF, attitude/orbit, etc.).
- **geometricIntegration**: Geometric integration methods (Lie group integrators, etc.).
- **gpr**: Gaussian Process Regression (GPR) functions.
- **hifiSRP**: High-fidelity Solar Radiation Pressure (SRP) modeling.
- **lightcurves**: Light curve calculation and analysis.
- **math**: General mathematical functions (Legendre polynomials, skew symmetric matrix, etc.).
- **object**: 3D object handling (reading .obj/.mtl files, self-shadowing).
- **orbit**: Orbital mechanics functions (Kepler's equation, coordinate transformations, EGM2008, EOP, etc.).
- **orbitDetermination**: Orbit determination algorithms (Gibbs, Gauss, Double R, etc.).
- **relativeOrbit**: Relative orbit dynamics and kinematics (HCW, ROE, etc.).
- **sphericalGaussian**: Spherical Gaussian functions.
- **srp**: Solar Radiation Pressure (SRP) models (Analytical, Cannonball, etc.).
- **sunMoon**: Sun and Moon position calculations, ephemeris.
- **time**: Time system conversions (JD, MJD, UTC, TT, etc.).
- **ukfCkf**: Unscented / Cubature Kalman Filter implementations.
- **utility**: General utility scripts (plotting, figure saving, etc.).

This repository uses Earth imagery from NASA's Blue Marble dataset.
Credit: NASA / Earth Observatory.