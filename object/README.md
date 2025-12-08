# abst.

衛星形状モデルとそれを任意の数の微小面に分割したものを扱う関数．BlenderやFusion等で3Dモデルを作成し，`obj`ファイルとして出力したものを読み込む，という想定をしています．

mtlファイルがある場合は，Ca, Cd, Csが自動的に読み込まれます．
次の場合は，`Ca = 0.960784, Cd = 0.960784, Cs = 0.964706`となります：
> newmtl アルミニウム_-_つや出し  
> Kd 0.960784 0.960784 0.964706

# note

`sat = readSC(’satFileName.obj’)` の構造体の構成（e.g., sat.normal, sat.faces）は次のとおり

微小面の数を`N`とする．単位は基本的に${\rm m, kg, kgm^{2}}$なので注意．`selfShadow.mlx`など，`sat`を使う関数をmex化するときはfieldの順番とsizeも大事になる（不用意に入れ替えたり，削除するとmex化できなくなる）．

| index | fields | size | 説明 | 単位 |
| --- | --- | --- | --- | --- |
| 1 | vertices | (> N) x 3 | 頂点の座標．なので 面の数より多くなる．<br>（e.g., 微小面が三角だとverticesが3つ繋がれて面が構成される．） | ${\rm m}$ |
| 2 | normal | N x 3 | 微小面の（外向き）法線ベクトル．単位ベクトル． | - |
| 3 | faces | N x 3 or N x 4 | 微小面を構成する`vertices`のindex．e.g., 1, 3, 5だと1番目，3番目，5番目のverticesを繋いで面を構成する． | - |
| 4 | area | N x 1 | 微小面の面積 | ${\rm m^{2}}$ |
| 5 | pos | N x 3 | 微小面の幾何中心位置ベクトル | ${\rm m}$ |
| 6 | uu | N x 3 | Ashikhmin–Shirleyモデル等を使うときの異方性基底ベクトル | - |
| 7 | uv | N x 3 | Ashikhmin–Shirleyモデル等を使うときの異方性基底ベクトル | - |
| 8 | qlb | N x 4 | 機体固定座標系から異方性基底ベクトルを$x,y$軸とするlocal frameへのquaternion | - |
| 9 | F0 | N x 1 | Fresnel reflection at normal incidence | - |
| 10 | kappa | N x 1 | 熱拡散係数 for SRP | - |
| 11 | Ca | N x 1 | 吸収率（mtlファイルがある場合は自動的に読み込まれる） | - |
| 12 | Cd | N x 1 | 拡散率（mtlファイルがある場合は自動的に読み込まれる） | - |
| 13 | Cs | N x 1 | 鏡面反射率（mtlファイルがある場合は自動的に読み込まれる） | - |
| 14 | nu | N x 1 | Ashikhmin–Shirleyモデルの異方性パラメータ | - |
| 15 | nv | N x 1 | Ashikhmin–Shirleyモデルの異方性パラメータ | - |
| 16 | mCT | N x 1 | Cook-Torranceモデルのroughnessパラメータ | - |
| 17 | fObs | N x 1 | for light curve | - |
| 18 | MOI | 3 x 3 | 慣性テンソル | ${\rm kgm^{2}}$ |
| 19 | m | 1 | 質量 | ${\rm kg}$ |
| 20 | sunlitFlag | N x 1 | self-shadowingのflag．0の時self-shadowingで影になっている | - |
| 21 | force | N x 3 | 機体固定座標系に対する外力．e.g., SRP | N |
| 22 | torque | N x 3 | 機体固定座標系に対する外力トルク．e.g., SRP torque | Nm |


# Functions

各関数の概要は以下の通りです．

- `readSC.m`
  `obj`ファイルを読み込み，構造体`sat`を作成するメイン関数です．
  面積や法線ベクトルなどの幾何学的情報を計算し，物理パラメータ（反射率など）も設定します．

- `selfShadow.m`
  衛星の各微小面に対して自己影（self-shadowing）の判定を行います．
  太陽方向ベクトルを入力とし，他の面によって遮蔽されている面を特定してフラグを立てます．

- `showSC.m`
  読み込んだ衛星モデル（`sat`構造体）を3次元プロットで表示します．
  各面の色分けや，影の可視化（`shadowFlag`の反映）など，形状確認に便利です．

- `calcRayIntersect.m`
  光線（Ray）と三角形（Triangle）の交差判定を行う計算関数です．
  `selfShadow`の内部で使用され，高速な交差判定アルゴリズムを実装しています．

- `calcAreaObj.m`
  各微小面の面積と幾何中心（重心）位置を計算します．
  `readSC`の内部で呼び出され，SRPや空力などの物理計算に必要な基本量を導出します．

- `calcLocalFrame.m`
  各面における局所座標系（facet上の基底ベクトルなど）を計算します．
  異方性反射モデル（Ashikhmin-Shirleyなど）を使用する場合に必要となります．

- `drawEarth.m`
  地球を3次元プロットで描画します．
  テクスチャ画像（`earth.jpg`）をマッピングし，GMSTに応じた自転も考慮して表示可能です．

# ref.

### 衛星モデル作成

[objファイル作成・読み込み](https://www.notion.so/obj-22b0a48416664f1081331654b3981a8a?pvs=21) 　を参照．

---

# English Description

# abst.

Functions for handling satellite shape models and their division into an arbitrary number of facets. It is assumed that 3D models are created using tools like Blender or Fusion, and loaded as exported `obj` files.

If an mtl file exists, Ca, Cd, and Cs are automatically loaded.
In the following case, `Ca = 0.960784, Cd = 0.960784, Cs = 0.964706`:
> newmtl Aluminum_-_Polished
> Kd 0.960784 0.960784 0.964706

# note

The structure configuration of `sat = readSC('satFileName.obj')` (e.g., sat.normal, sat.faces) is as follows:

Let `N` be the number of facets. Note that units are basically ${\rm m, kg, kgm^{2}}$. When mexing functions that use `sat`, such as `selfShadow.mlx`, the order and size of fields are important (randomly swapping or deleting them will prevent mexing).

| index | fields | size | Description | Unit |
| --- | --- | --- | --- | --- |
| 1 | vertices | (> N) x 3 | Vertex coordinates. Thus, more than the number of faces. <br>(e.g., if facets are triangles, 3 vertices are connected to form a face.) | ${\rm m}$ |
| 2 | normal | N x 3 | Normal vector of the facet (outward). Unit vector. | - |
| 3 | faces | N x 3 or N x 4 | Indices of `vertices` forming the facet. e.g., 1, 3, 5 means connecting the 1st, 3rd, and 5th vertices to form a face. | - |
| 4 | area | N x 1 | Area of the facet | ${\rm m^{2}}$ |
| 5 | pos | N x 3 | Geometric center position vector of the facet | ${\rm m}$ |
| 6 | uu | N x 3 | Anisotropic basis vector when using Ashikhmin–Shirley model etc. | - |
| 7 | uv | N x 3 | Anisotropic basis vector when using Ashikhmin–Shirley model etc. | - |
| 8 | qlb | N x 4 | Quaternion from body-fixed frame to local frame with anisotropic basis vectors as x, y axes | - |
| 9 | F0 | N x 1 | Fresnel reflection at normal incidence | - |
| 10 | kappa | N x 1 | Thermal diffusion coefficient for SRP | - |
| 11 | Ca | N x 1 | Absorption coefficient (automatically loaded if mtl file exists) | - |
| 12 | Cd | N x 1 | Diffuse reflection coefficient (automatically loaded if mtl file exists) | - |
| 13 | Cs | N x 1 | Specular reflection coefficient (automatically loaded if mtl file exists) | - |
| 14 | nu | N x 1 | Anisotropic parameter for Ashikhmin–Shirley model | - |
| 15 | nv | N x 1 | Anisotropic parameter for Ashikhmin–Shirley model | - |
| 16 | mCT | N x 1 | Roughness parameter for Cook-Torrance model | - |
| 17 | fObs | N x 1 | for light curve | - |
| 18 | MOI | 3 x 3 | Inertia tensor | ${\rm kgm^{2}}$ |
| 19 | m | 1 | Mass | ${\rm kg}$ |
| 20 | sunlitFlag | N x 1 | Flag for self-shadowing. 0 means shadowed by self-shadowing. | - |
| 21 | force | N x 3 | External force in body-fixed frame. e.g., SRP | N |
| 22 | torque | N x 3 | External torque in body-fixed frame. e.g., SRP torque | Nm |


# Functions

The overview of each function is as follows:

- `readSC.m`
  Main function to load `obj` files and create the `sat` structure.
  It calculates geometric information such as area and normal vectors, and also sets physical parameters (reflectivity, etc.).

- `selfShadow.m`
  Determines self-shadowing for each facet of the satellite.
  Takes the sun direction vector as input, identifies faces occluded by others, and sets flags.

- `showSC.m`
  Displays the loaded satellite model (`sat` structure) in a 3D plot.
  Useful for shape verification, including color-coding of faces and shadow visualization (reflecting `shadowFlag`).

- `calcRayIntersect.m`
  Calculation function for ray-triangle intersection testing.
  Used internally by `selfShadow`, implementing a fast intersection algorithm.

- `calcAreaObj.m`
  Calculates the area and geometric center (centroid) of each facet.
  Called internally by `readSC` to derive basic quantities required for physical calculations such as SRP and aerodynamics.

- `calcLocalFrame.m`
  Calculates the local coordinate system (basis vectors on the facet, etc.) for each face.
  Required when using anisotropic reflection models (such as Ashikhmin-Shirley).

- `drawEarth.m`
  Draws the Earth in a 3D plot.
  Maps a texture image (`earth.jpg`) and can display rotation according to GMST.

# ref.

### Satellite Model Creation

Refer to [Creating and Loading obj files](https://www.notion.so/obj-22b0a48416664f1081331654b3981a8a?pvs=21).
