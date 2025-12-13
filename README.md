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
