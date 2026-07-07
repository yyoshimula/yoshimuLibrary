# yoshimuLibrary マニュアル

吉村が宇宙機の姿勢・軌道力学研究の中で作成・整備している MATLAB ライブラリ **yoshimuLibrary** のリファレンスマニュアルです。各関数のソースコード（`%[text]` ドキュメントコメントと `function` シグネチャ）を突き合わせて記述しており、末尾の「[付録: コードとドキュメントの不整合](#付録-コードとドキュメントの不整合-要確認)」に、コードと元コメントの食い違いをまとめています。

> バグ・記述ミスを見つけたら → y.yoshimura.a64@m.kyushu-u.ac.jp / y.yoshimula@gmail.com

## 使い方

一時的に使う場合は対象ディレクトリに `addpath` でパスを通します。継続的に使う場合は MATLAB の「パスの設定 (Set Path)」で本リポジトリ配下を追加してください。

```matlab
% 例: 軌道・姿勢・SRP を使う
addpath(genpath('/path/to/yoshimuLibrary'));   % 一時的に再帰追加
```

リポジトリの取得（Git LFS 必須 — 係数テーブルやテクスチャが LFS 管理）:

```bash
git clone git@github.com:yyoshimula/yoshimuLibrary.git
cd yoshimuLibrary
git lfs install
git lfs pull
```

- **動作要件**: MATLAB R2024a 以降、Git LFS。
- 一部関数は係数ファイル（EGM2008 / IAU2006 / ELP / VSOP / JB2008 / スペースウェザー等）を同梱データから読み込みます。該当ディレクトリごと addpath しておくこと。

## 規約（全関数共通の約束事）

- **引数は原則「横ベクトル (row vector)」**。ODE で運動方程式を解くと「行方向＝時系列、列方向＝状態変数」の並びになるため、これに合わせています。多くの関数は `n×3` などの行列を受け取り、**時系列データを for 文なしで一括変換**できるようベクトル化されています（`n` が時刻サンプル数）。
- **ドキュメントはコード内コメントが正**。各 `.m` は先頭に `%[text]` 形式のコメント（タイトル・inputs・outputs・note・references・See also）を持ち、MATLAB 上では Live Script として整形表示できます。本マニュアルはこのコメントと実シグネチャの両方から生成しています。
- **`.m` が保守対象**。かつては同一処理の Live 関数 `****.mlx` と m ファイル `****_.m` を併存させていましたが、現在リポジトリはほぼ `.m` に一本化されています（`.mlx` は2ファイルのみ）。引数検証の `arguments` ブロックは多くがコメントアウトされており、計算速度を優先した状態です。
- **依存する主な構造体**:
  - `sat` — 宇宙機形状モデル。`object/readSC` で `.obj`/`.mtl` から読み込む。面法線・面積・BRDF パラメータ等を保持し、`srp` / `lightcurves` / `hifiSRP` が入力に取る。
  - `const` — 定数構造体。`orbit/orbitConst` で生成（地球重力定数・半径・光速・AU 等）。
  - `oe` — 軌道要素 `[a, e, i, Ω, ω, f または M]`（`orbit` / `relativeOrbit` 系）。
- **角度は原則 rad、距離の単位は関数内で統一必須**（重力定数・位置・速度・軌道長半径の単位を揃えること）。

---

## 目次


1. [数学ユーティリティ (`math/`)](#数学ユーティリティ-math)
2. [単位・角度変換 (`conversion/`)](#単位角度変換-conversion)
3. [時刻系変換 (`time/`)](#時刻系変換-time)
4. [軌道力学・座標変換 (`orbit/`)](#軌道力学座標変換-orbit)
5. [初期軌道決定 (`orbitDetermination/`)](#初期軌道決定-orbitdetermination)
6. [相対軌道 (ROE / 視線角) (`relativeOrbit/`)](#相対軌道-roe--視線角-relativeorbit)
7. [太陽・月の位置と重力 (`sunMoon/`)](#太陽月の位置と重力-sunmoon)
8. [宇宙環境モデル (IGRF・大気密度) (`environment/`)](#宇宙環境モデル-igrf大気密度-environment)
9. [姿勢表現・運動学 (`attitude/`)](#姿勢表現運動学-attitude)
10. [デュアルクォータニオン (`dualQuaternions/`)](#デュアルクォータニオン-dualquaternions)
11. [幾何学的積分 (`geometricIntegration/`)](#幾何学的積分-geometricintegration)
12. [宇宙機形状モデル (OBJ) (`object/`)](#宇宙機形状モデル-obj-object)
13. [太陽輻射圧 (SRP) (`srp/`)](#太陽輻射圧-srp-srp)
14. [高精度 SRP 近似 (`hifiSRP/`)](#高精度-srp-近似-hifisrp)
15. [ライトカーブ (`lightcurves/`)](#ライトカーブ-lightcurves)
16. [球面ガウス関数 (`sphericalGaussian/`)](#球面ガウス関数-sphericalgaussian)
17. [UKF / CKF フィルタ (`ukfCkf/`)](#ukf--ckf-フィルタ-ukfckf)
18. [ガウス過程回帰 (GPR) (`gpr/`)](#ガウス過程回帰-gpr-gpr)
19. [汎用ユーティリティ (`utility/`)](#汎用ユーティリティ-utility)
20. [可視化 (`visualization/`)](#可視化-visualization)
21. [サンプルスクリプト (`examples/`)](#サンプルスクリプト-examples)
22. [補助ディレクトリ（SPICE・演習・テスト・参考資料）](#補助ディレクトリspice演習テスト参考資料)
23. [付録: コードとドキュメントの不整合 (要確認)](#付録-コードとドキュメントの不整合-要確認)


---

## 数学ユーティリティ (`math/`)

姿勢・軌道計算で共通して使う基本的な数学プリミティブを集めたディレクトリです。ベクトルの外積行列化、角度の正規化、行方向に並んだベクトル群の正規化、球面調和関数で用いるルジャンドル陪多項式を提供します。`normRow` は行方向 = 1 サンプル・列方向 = 成分（横ベクトル）を各行とする行列を想定し、各行を独立に単位化します。外部の構造体（`sat` / `const` / `oe` 等）への依存はありません。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `skew` | 3 次元ベクトルから外積（歪対称）行列を生成 |
| `wrapPi` | 角度をラジアンで [-π, π) に正規化 |
| `normRow` | 行列の各行ベクトルを 2-ノルムで正規化 |
| `associatedLegendre` | ルジャンドル陪多項式 P_l^m(x) を全位数 m=0..l で計算 |

### Core

#### `skew`
```matlab
S = skew(x)
```
3 次元ベクトル x から外積行列（歪対称行列）を生成します。任意のベクトル a に対して `skew(x)*a` が `cross(x, a)` に等しくなります。
- **入力**: `x` (3×1 または 1×3) — 3 次元ベクトル
- **出力**: `S` (3×3) — 歪対称行列 `[0 -x3 x2; x3 0 -x1; -x2 x1 0]`

#### `wrapPi`
```matlab
lambda = wrapPi(lambda)
```
角度（ラジアン）を半開区間 [-π, π) に折り返します。配列入力は要素ごとに処理されます。
- **入力**: `lambda` (任意サイズ, rad) — 折り返す角度（配列可）
- **出力**: `lambda` (入力と同サイズ, rad) — [-π, π) に正規化された角度
- **Note**: 出力区間は半開 [-π, π)。実装は `mod(lambda + pi, 2*pi) - pi` で、+π ちょうどは -π に写る。
- *See also*: `mod`

#### `normRow`
```matlab
B = normRow(A)
```
行列 A の各行を 2-ノルムで割って単位ベクトル化します（行方向 = サンプル、列方向 = 成分の横ベクトル規約）。
- **入力**: `A` (N×D) — 行ごとに正規化したいベクトルを並べた行列
- **出力**: `B` (N×D) — 各行を単位化した行列（ノルムが eps 以下の場合は入力をそのまま返す）
- **Note**: ゼロ割回避の分岐 `if (Anorm > eps)` は、`Anorm` が行ごとのノルムからなる列ベクトル（N>1）のとき **全行が eps 超のときだけ** 正規化ブランチに入る。ゼロ行と非ゼロ行が混在すると条件全体が偽になり、行列全体が未正規化のまま返る。単一行または全行非ゼロでの使用を想定。

### 球面調和・特殊関数

#### `associatedLegendre`
```matlab
[P, m_values] = associatedLegendre(l, x)
```
次数 l のルジャンドル陪多項式 P_l^m(x) を、全位数 m = 0, 1, ..., l についてまとめて計算します。内部補助関数 `legendreRecursive`（漸化式）と `double_factorial`（二重階乗）を用います。
- **入力**: `l` (スカラー) — 次数（非負整数）; `x` (スカラーまたはベクトル) — 評価点
- **出力**: `P` ((l+1)×length(x)) — `P(i,j) = P_l^{m_i}(x_j)`（第 i 行が位数 m=i-1）; `m_values` (1×(l+1)) — 位数ベクトル [0, 1, ..., l]
- **Note**: 実装は Condon–Shortley 位相 `(-1)^m` を含む **非正規化**の標準漸化式。ファイル冒頭コメントは「正規化された」「複数の方法を実装」と書くが、実際には正規化は行われず漸化法 1 種のみ（下の不整合参照）。`legendreRecursive` 内の入力チェック（m の範囲・|x|>1 警告）はコメントアウトされており無効。

---

## 単位・角度変換 (`conversion/`)

物理量の単位変換と角度表現の相互変換を行う小関数群。すべて要素ごと演算（`.*` `./`）なのでスカラでも配列でも同じ形状で返り、時系列・状態ベクトルのどの向きに並べても使える。距離系の AU⇔km 変換は天文単位の定義値を保持する `const` 構造体（`const.AU`、`orbit/orbitConst` が生成）に依存する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `arcs2rad` | 角度を秒角 (arcsec) からラジアンへ変換 |
| `rad2arcs` | 角度をラジアンから秒角 (arcsec) へ変換 |
| `hms2deg` | 時・分・秒で表した角度を度 (deg) へ変換 |
| `au2km` | 距離を天文単位 (AU) から km へ変換 |
| `km2au` | 距離を km から天文単位 (AU) へ変換 |
| `s2day` | 時間を秒 (s) から日 (day) へ変換 |

### 角度変換

#### `arcs2rad`
```matlab
rad = arcs2rad(angle)
```
秒角 (arcsec) をラジアンに変換する。
- **入力**: `angle` (任意サイズ, arcsec) — 変換する角度
- **出力**: `rad` (`angle` と同サイズ, rad) — 変換後の角度

#### `rad2arcs`
```matlab
arcs = rad2arcs(angle)
```
ラジアンを秒角 (arcsec) に変換する。`arcs2rad` の逆変換。
- **入力**: `angle` (任意サイズ, rad) — 変換する角度
- **出力**: `arcs` (`angle` と同サイズ, arcsec) — 変換後の角度
- **Note**: ヘッダコメントのタイトルおよび出力の記述が `arcs2rad` からの流用で誤っている（下記「不整合」参照）。実装は rad→arcsec で正しい。
- *See also*: `arcs2rad`

#### `hms2deg`
```matlab
out = hms2deg(hour, min, sec)
```
時・分・秒で表した角度を度 (deg) に換算する。
- **入力**: `hour` (任意サイズ, hour) — 時成分; `min` (同, min) — 分成分; `sec` (同, sec) — 秒成分
- **出力**: `out` (入力と同サイズ, deg) — 変換後の角度
- **Note**: 実装は `out = hour*15 + (min*60 + sec)/3600`。標準的な時角→度の換算（`15*(hour + min/60 + sec/3600)`）とは分・秒項の係数が一致しないため、分・秒成分を含む入力では想定と異なる値になる可能性がある（要確認）。

### 距離変換

#### `au2km`
```matlab
km = au2km(AU, const)
```
距離を天文単位 (AU) から km に変換する。
- **入力**: `AU` (任意サイズ, AU) — 変換する距離; `const` (構造体) — `const.AU`（1 AU の km 値）を持つ定数構造体
- **出力**: `km` (`AU` と同サイズ, km) — 変換後の距離
- **Note**: `const.AU` は `orbit/orbitConst` が生成する定数を渡す前提。
- *See also*: `km2au`, `orbitConst`

#### `km2au`
```matlab
AU = km2au(km, const)
```
距離を km から天文単位 (AU) に変換する。`au2km` の逆変換。
- **入力**: `km` (任意サイズ, km) — 変換する距離; `const` (構造体) — `const.AU` を持つ定数構造体
- **出力**: `AU` (`km` と同サイズ, AU) — 変換後の距離
- **Note**: `const.AU` は `orbit/orbitConst` が生成する定数を渡す前提。
- *See also*: `au2km`, `orbitConst`

### 時間変換

#### `s2day`
```matlab
day = s2day(s)
```
時間を秒 (s) から日 (day) に変換する（`s/86400`）。
- **入力**: `s` (任意サイズ, s) — 変換する時間
- **出力**: `day` (`s` と同サイズ, day) — 変換後の時間

---

## 時刻系変換 (`time/`)

暦日・ユリウス日・修正ユリウス日の相互変換と、UTC/TAI/TT/UT1 といった時刻系間のオフセット（うるう秒 ΔAT や ΔT）を扱う関数群。多くの関数は入力を列ベクトル `(:,1)` として受け取り、複数エポックをまとめて処理できる。うるう秒処理は `leapS`（テーブル生成）→ `dAT`（各 JD の ΔAT 計算）の2段構成で、`utc2tt` はその ΔAT を受け取って TT を返す。秒⇔日の換算には `conversion/` の `s2day` に依存する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `gc2jd` | グレゴリオ暦日をユリウス日 (JD) に変換 |
| `jd2gc` | ユリウス日 (JD) をグレゴリオ暦日に変換 |
| `doy2gc` | 通日 (day of year) をグレゴリオ暦日に変換 |
| `jd2mjd` | ユリウス日 (JD) を修正ユリウス日 (MJD) に変換 |
| `mjd2jd` | 修正ユリウス日 (MJD) をユリウス日 (JD) に変換 |
| `jd2jdT` | ユリウス日を J2000.0 からのユリウス世紀数 T に変換 |
| `leapS` | うるう秒テーブル（有効開始 JD と累積うるう秒）を生成 |
| `dAT` | 各 JD の ΔAT (= TAI − UTC) を計算 |
| `utc2tt` | UTC (JD) と ΔAT から TT (JD) を計算 |
| `ut2tt` | UT1 から TT へのオフセット Δt (= TT − UT1) を計算 |

### 暦日 ⇔ ユリウス日

#### `gc2jd`
```matlab
jd = gc2jd(year, month, day, hour, minute, second)
```
グレゴリオ暦日をユリウス日 (JD) に変換する。各引数は同じ長さの列ベクトルにでき、まとめて変換できる。
- **入力**: `year` (n×1) — 年; `month` (n×1) — 月; `day` (n×1) — 日; `hour` (n×1) — 時; `minute` (n×1) — 分; `second` (n×1) — 秒
- **出力**: `jd` (n×1, day) — ユリウス日
- **Note**: Meeus, *Astronomical Algorithms* (1998) p.61 Eq.(7.1) に基づく。
- *See also*: `jd2gc`

#### `jd2gc`
```matlab
[year, month, day, hour, minute, second] = jd2gc(jd)
```
ユリウス日 (JD) からグレゴリオ暦日（年月日時分秒）を復元する。
- **入力**: `jd` (n×1, day) — ユリウス日
- **出力**: `year`, `month`, `day`, `hour`, `minute`, `second`（各 n×1）— 暦日成分

#### `doy2gc`
```matlab
[month, day, hour, minute, second] = doy2gc(year, doy)
```
通日（year の 1 月 1 日を 1.0 とする小数の day of year）をグレゴリオ暦日に変換する。内部で `gc2jd` と `jd2gc` を経由する。
- **入力**: `year` — 年; `doy` — 通日（1 始まり、1.0 〜 366.x の小数）
- **出力**: `month`, `day`, `hour`, `minute`, `second` — 暦日成分（年は入力と同じなので返さない）
- *See also*: `gc2jd`, `jd2gc`

### ユリウス日系の変換

#### `jd2mjd`
```matlab
mjd = jd2mjd(jd)
```
ユリウス日 (JD) を修正ユリウス日 (MJD = JD − 2400000.5) に変換する。
- **入力**: `jd` (n×1, day) — ユリウス日
- **出力**: `mjd` (n×1, day) — 修正ユリウス日

#### `mjd2jd`
```matlab
jd = mjd2jd(mjd)
```
修正ユリウス日 (MJD) をユリウス日 (JD = MJD + 2400000.5) に変換する。
- **入力**: `mjd` (n×1, day) — 修正ユリウス日
- **出力**: `jd` (n×1, day) — ユリウス日

#### `jd2jdT`
```matlab
T = jd2jdT(jd)
```
ユリウス日を J2000.0 (JD 2451545.0) を基準としたユリウス世紀数 T = (jd − 2451545.0) / 36525 に変換する。
- **入力**: `jd` (n×1, day) — ユリウス日
- **出力**: `T` (n×1) — ユリウス世紀数（36525 日 = 1 ユリウス世紀）

### うるう秒・時刻系オフセット

#### `leapS`
```matlab
leapJD = leapS()
```
MATLAB 組み込みの `leapseconds` テーブルから、うるう秒の有効開始 JD と累積うるう秒の対応表を生成する。
- **入力**: なし
- **出力**: `leapJD` (m×2) — 1 列目が各うるう秒の有効開始 JD、2 列目がその時点での累積うるう秒（TAI − UTC の増分累積）
- **Note**: うるう秒は挿入日末尾（23:59:60 UTC）に入るため、新しい TAI − UTC は翌日 00:00 UTC から有効。`leapJD(:,1)` はこの有効開始日の JD（挿入日 +1 日）。結果は `dAT` の引数 `leapJD` に渡して使う。
- *See also*: `dAT`

#### `dAT`
```matlab
deltaAT = dAT(jd, leapJD)
```
各エポックにおける ΔAT (= TAI − UTC) を計算する。`leapJD` 表の有効開始 JD が対象 JD 以前に何件あるかを数え、基準値 10 秒（1972-01-01 以降）に累積うるう秒を加える。
- **入力**: `jd` (スカラまたはベクトル, day) — ユリウス日; `leapJD` (m×2) — `leapS` の出力（有効開始 JD と累積うるう秒）
- **出力**: `deltaAT` (jd と同サイズ, s) — ΔAT = TAI − UTC
- **Note**: 使用前に `leapS` を実行し、その結果を `leapJD` に渡すこと。`deltaAT` は `jd` と同じ形状で返る。Vallado & McClain (2001) に基づく。
- *See also*: `leapS`

#### `utc2tt`
```matlab
jdTT = utc2tt(jdUTC, deltaAT)
```
UTC のユリウス日と ΔAT から TT (Terrestrial Time) のユリウス日を計算する（TT = UTC + ΔAT + 32.184 s）。
- **入力**: `jdUTC` (n×1, day) — UTC のユリウス日; `deltaAT` (n×1, s) — ΔAT (= TAI − UTC)、`dAT` で得る値
- **出力**: `jdTT` (n×1, day) — TT のユリウス日
- **Note**: 秒→日の換算に `conversion/s2day` を使う。ΔAT は `dAT` で別途求めて渡す。Vallado & McClain (2001) p.220。

#### `ut2tt`
```matlab
dt = ut2tt(jd)
```
UT1 から TT へのオフセット Δt (= TT − UT1) を、年をキーとした経験的モデル（1620 年以降の実測補間 + 前後の多項式外挿）で計算する。
- **入力**: `jd` (n×1, day) — ユリウス日
- **出力**: `dt` (n×1, s) — Δt = TT − UT1
- **Note**: NASA TP-2006-214141 (Espenak & Meeus) の ΔT 近似式。1620.5〜2019.5 年は3年ごとの実測テーブルを2次補間、2019.5〜2050 年・2150 年未満・それ以降でそれぞれ別の多項式に切り替わるため、区間境界で連続性が完全には保証されない。ΔAT のようなうるう秒の整数ステップではなく地球回転の変動を含む量である点に注意。

---

## 軌道力学・座標変換 (`orbit/`)

軌道要素・状態ベクトルの相互変換、Kepler 方程式・各種近点角、Gauss/変分方程式による摂動、地球重力（EGM2008）と地球影、そして地球回転（歳差・章動・極運動）を含む慣性系↔地球固定系の座標・フレーム変換をまとめたモジュール。多くの関数は時系列データを行方向に並べる横ベクトル規約（例: 位置 `nx3`、Julian day `nx1`）を採る。共通の定数は構造体 `const`（`orbitConst` が生成）で受け渡し、EGM2008 係数は構造体 `EGM`、VSOP87 係数は `earthVSOP` にまとめる。軌道要素をフィールド展開した構造体 `oe`（`.a/.e/.inc/.raan/.w/.nu/...`）を入出力に使う関数群（`calcOrbitalState`, `coe2mee`, `mee2coe`, `gve`, `mee`）と、単純なベクトル `[a,e,i,Ω,w,ν]` を使う関数群（`oe2rv`, `rv2oe`）が混在する点に注意。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `oe2rv` | 軌道要素 `[a,e,i,Ω,w,ν(またはM)]` から位置・速度ベクトルへ変換 |
| `rv2oe` | 位置・速度ベクトルから軌道要素 `[a,e,i,Ω,w,ν]` へ変換 |
| `calcOrbitalState` | `oe` 構造体に p/h/n/u/M/r と位置・速度ベクトルを補完 |
| `coe2mee` | 古典軌道要素 → 修正春分点要素 (MEE) |
| `mee2coe` | 修正春分点要素 (MEE) → 古典軌道要素 |
| `dcmI2RTN` | 慣性系から RTN 系への方向余弦行列 |
| `gve` | Gauss 変分方程式による古典軌道要素の時間微分 |
| `mee` | 変分方程式による修正春分点要素の時間微分 |
| `eAnomaly` | 離心率・真近点角から離心近点角 E |
| `meanAnomaly` | 離心率・真近点角から平均近点角 M |
| `trueAnomaly` | a,e,M から真近点角 f と動径距離 r |
| `keplerEq` | Kepler 方程式を Newton 法で解き離心近点角 E |
| `mean2Osc` | 平均軌道要素 → 接触軌道要素 |
| `earthG` | ECEF 系での地球重力加速度（EGM2008、SPICE 併用可） |
| `egm2008` | EGM2008 ジオポテンシャルによる非球対称地球重力（ECEF, デカルト） |
| `readEGM2008` | EGM2008 係数ファイルを読み `EGM.Cnm/.Snm` を返す |
| `shadow` | 地球影関数（1: 日照、0: 本影） |
| `earthVSOP87` | VSOP87 による地球の日心黄経・黄緯・距離 |
| `vsopConst` | VSOP87 の係数テーブル `earthVSOP` (スクリプト的定数関数) |
| `orbitConst` | 軌道伝播用の定数構造体 `const` |
| `ecef2LatLonH` | ECEF 位置から測地緯度・経度・高度 |
| `geodetic2Geocentric` | 測地緯度・経度・高度から ECEF 位置 |
| `geocentric2Geodetic` | 地心直交座標から測地緯度・経度・高度 |
| `gmst` | Julian day から Greenwich 平均恒星時 |
| `gast` | Julian day から Greenwich 視恒星時 (IAU-76/FK5) |
| `era` | 地球回転角 (IAU-2006/2000, CIO ベース) |
| `earthW` | LOD から地球自転角速度 |
| `wobble` | 極運動の年周ゆらぎ s'（IAU-2006/2000, CIO ベース） |
| `obliquity` | 与えられた元期の黄道傾斜角（平均） |
| `nutation` | 平均→真赤道 (TOD) の章動量 dψ, dε |
| `nutationDCM` | MOD→TOD の章動 DCM (IAU-76/FK5) |
| `nutationQ` | MOD→TOD の章動クォータニオン |
| `precession` | 平均→TOD の歳差角一式 |
| `precessionDCM` | J2000(ICRF)→MOD の歳差 DCM |
| `precessionQ` | 歳差クォータニオン |
| `precessionNutation` | IAU-2006/2000 理論による歳差＋章動 DCM (CIRS→GCRF) |
| `earthFullRotQ` | 歳差＋章動を含む地球回転クォータニオン |
| `earthNutationPrecessionQ` | 歳差＋章動を含む地球回転クォータニオン（極運動なし） |
| `itrf2gcrf` | ITRF→GCRF 回転行列 (IAU-2006/2000, CIO) |
| `qITRF2gcrf` | ITRF→GCRF クォータニオン（`itrf2gcrf` のラッパ） |
| `pef2itrf` | PEF→ITRF の DCM（極運動） |
| `mod2J2000` | MOD 系の軌道要素を J2000 系へ変換 |
| `teme2J2000` | TEME 系の軌道要素を J2000(FK5) 系へ変換 |
| `tod2Mod` | TOD 系の軌道要素を MOD 系へ変換 |
| `eop` | 指定日付の地球姿勢パラメータ (EOP) を抽出 |
| `readEOP` | EOP データファイルを読み込み |
| `readIAU06` | IAU06 (x,y,z) 係数を読み込み |
| `readJB2008` | Jacchia–Bowman 2008 の係数群を読み込み |
| `readTLE` | TLE ファイルを読み平均軌道要素等を返す |
| `compareLegendre` | 随伴 Legendre 関数の 2 実装を比較 (スクリプト) |

### 軌道要素と状態ベクトル

#### `oe2rv`
```matlab
[r, v] = oe2rv(oe, flag, mu)
```
軌道要素から慣性系の位置・速度ベクトルを計算する。
- **入力**: `oe` (nx6) — `[a, e, i, Ω, w, f(またはM)]`（長さ: m または km、-、rad×4）; `flag` (1,1) — 1: 第6要素は真近点角、0: 平均近点角; `mu` (1,1) — 中心天体の重力定数（`r,v,a` と単位を統一）
- **出力**: `r` (nx3, m または km) — 位置; `v` (nx3, m/s または km/s) — 速度
- **Note**: 重力定数・位置速度・軌道長半径は単位を揃えること。円軌道 (e≈0) では w を、赤道軌道 (i≈0 または π) では Ω を 0 に丸めて特異を回避する。

#### `rv2oe`
```matlab
oe = rv2oe(r, v, mu)
```
位置・速度ベクトルから古典軌道要素を計算する。円/赤道軌道は Vallado の rv2coe 規約に従う。
- **入力**: `r` (nx3, km) — 位置; `v` (nx3, km/s) — 速度; `mu` (1,1) — 地球重力定数
- **出力**: `oe` (nx6) — `[a, e, i, Ω, w, ν]`（km、-、rad×4）
- **Note**: 出力は 6 列の行列であり、doc コメントが列挙する M・u・trueLon は内部計算のみでコード上は返さない。円軌道では w=0 とし ν を緯度引数（円赤道では真経度）に置換、赤道軌道では Ω=0 とする。放物線軌道 (|ξ|≤eps) は a=Inf。

#### `calcOrbitalState`
```matlab
oe = calcOrbitalState(oe, mu)
```
軌道要素構造体 `oe` に派生量（半通径 p、角運動量 h、平均運動 n、緯度引数 u、平均近点角 M、動径 r 等）と位置・速度ベクトルを補完して返す。
- **入力**: `oe` (struct) — 少なくとも `.a/.e/.inc/.raan/.w/.nu` を持つ構造体; `mu` (1,1) — 重力定数
- **出力**: `oe` (struct) — 上記に `.p/.h/.n/.u/.M/.uM/.r/.rVec/.vVec` を追加した構造体
- **Note**: `gve` から前処理として呼ばれる。内部で `oe2rv`, `meanAnomaly` を利用。

#### `coe2mee`
```matlab
oe = coe2mee(oe)
```
古典軌道要素から修正春分点要素 (MEE) `[p, f, g, h, k, L]` を計算し、`oe` 構造体のフィールド `.p_/.f_/.g_/.h_/.k_/.L_` に格納する。
- **入力**: `oe` (struct) — `.a/.e/.inc/.raan/.w/.nu` を持つ構造体
- **出力**: `oe` (struct) — 上記に MEE フィールド `.p_/.f_/.g_/.h_/.k_/.L_`（rad）を追加

#### `mee2coe`
```matlab
oe = mee2coe(oe)
```
修正春分点要素から古典軌道要素を計算し、`oe` 構造体のフィールド `.a/.e/.inc/.raan/.w/.nu` に格納する。
- **入力**: `oe` (struct) — MEE フィールド `.p_/.f_/.g_/.h_/.k_/.L_` を持つ構造体
- **出力**: `oe` (struct) — 古典軌道要素フィールドを追加

#### `dcmI2RTN`
```matlab
R = dcmI2RTN(raan, inc, w, nu)
```
慣性系から RTN（動径・接線・法線）系への方向余弦行列を Z-X-Z 回転（Ω, i, w+ν）で構成する。
- **入力**: `raan` (rad) — 昇交点赤経; `inc` (rad) — 軌道傾斜角; `w` (rad) — 近点引数; `nu` (rad) — 真近点角
- **出力**: `R` (3x3) — 慣性系→RTN 系の DCM

### 摂動・変分方程式

#### `gve`
```matlab
dOEdt = gve(oe, aRTN, anomalyFlag, mu)
```
Gauss 変分方程式で古典軌道要素の時間微分を計算する。内部で `calcOrbitalState` を呼び派生量を補完する。
- **入力**: `oe` (struct) — 軌道要素構造体; `aRTN` (3要素) — RTN 系の摂動加速度 [a_R, a_T, a_N]; `anomalyFlag` (1,1) — 1: 真近点角基準、それ以外: 平均近点角基準; `mu` (1,1) — 重力定数
- **出力**: `dOEdt` (6x1) — `[da/dt, de/dt, di/dt, dΩ/dt, dw/dt, dν(またはM)/dt]`

#### `mee`
```matlab
dOEdt = mee(oe, aRTN, mu)
```
変分方程式で修正春分点要素の時間微分を計算する。
- **入力**: `oe` (struct) — MEE フィールド `.p_/.f_/.g_/.h_/.k_/.L_` を持つ構造体; `aRTN` (3要素) — RTN 系摂動加速度 [a_R, a_T, a_N]; `mu` (1,1) — 重力定数
- **出力**: `dOEdt` (6x1) — `[dp/dt, df/dt, dg/dt, dh/dt, dk/dt, dL/dt]`

### 近点角・Kepler 方程式

#### `eAnomaly`
```matlab
E = eAnomaly(e, f)
```
離心率と真近点角から離心近点角を計算する。
- **入力**: `e` — 離心率; `f` (rad) — 真近点角
- **出力**: `E` (rad) — 離心近点角

#### `meanAnomaly`
```matlab
M = meanAnomaly(e, f)
```
離心率と真近点角から平均近点角を計算する。
- **入力**: `e` — 離心率; `f` (rad) — 真近点角
- **出力**: `M` (rad) — 平均近点角

#### `trueAnomaly`
```matlab
[f, r] = trueAnomaly(a, e, M)
```
平均近点角から真近点角と地心距離を計算する（内部で Kepler 方程式を解く）。
- **入力**: `a` (nx1, km) — 軌道長半径; `e` (nx1) — 離心率; `M` (nx1, rad) — 平均近点角
- **出力**: `f` (nx1, rad) — 真近点角; `r` (nx1, km) — 地心距離

#### `keplerEq`
```matlab
E = keplerEq(M, e, TOL)
```
Kepler 方程式を Newton 法で解き離心近点角を求める。
- **入力**: `M` (nx1, rad) — 平均近点角; `e` (nx1) — 離心率; `TOL` — 収束許容誤差（省略時 1e-8）
- **出力**: `E` (rad) — 離心近点角

#### `mean2Osc`
```matlab
osc = mean2Osc(n, e, i, Ome, w, M, const)
```
平均軌道要素を接触軌道要素へ変換する（Vallado の J2 短周期補正）。
- **入力**: `n` (rev/day) — 平均運動; `e` — 離心率; `i` (rad) — 軌道傾斜角; `Ome` (rad) — 昇交点経度; `w` (rad) — 近点引数; `M` (rad) — 平均近点角; `const` (struct) — `orbitConst` 定数
- **出力**: `osc` — 接触軌道要素ベクトル `[aOsc, eOsc, iOsc, OmeOsc, wOsc, fOsc, MOsc, rOsc, drOsc, pOsc, uOsc]`
- **Note**: doc の入力欄には `f` の記載もあるが、シグネチャは平均近点角 `M` を受け取る（第6引数）。

### 地球重力・地球影

#### `earthG`
```matlab
aEarth = earthG(jd, rVec, const, EGM, options)
```
ECEF 系での地球重力加速度を EGM2008 で計算する。既定では自前の歳差・章動・GAST 回転、`options.SPICE='on'` で SPICE の座標変換を使う。
- **入力**: `jd` (nx1, day) — Julian day; `rVec` (nx3, km) — 慣性系での衛星位置; `const` (struct) — 軌道定数; `EGM` (struct) — 重力係数（`.GEODEG/.Cnm/.Snm`）; `options.SPICE` — `'on'`/`'off'`（既定 `'off'`）
- **出力**: `aEarth` (nx3, km/s^2) — ECEF 系の地球重力加速度
- **Note**: doc の出力欄は「ECEF 系」と記すが、内部では ECEF で計算後に慣性系へ戻して返す実装（利用時に確認推奨）。

#### `egm2008`
```matlab
a = egm2008(rVec, deg, Cnm, Snm, const)
```
EGM2008 ジオポテンシャル（最大 20 次数まで）による非球対称地球重力の摂動加速度を ECEF・デカルト座標で計算する。
- **入力**: `rVec` (1x3, km) — ECEF 系での位置; `deg` — 使用する次数・位数; `Cnm`, `Snm` — EGM2008 の非正規化係数; `const` (struct) — 軌道定数
- **出力**: `a` (1x3, km/s^2) — ECEF 系（デカルト）の摂動加速度
- **Note**: 非正規化係数を使う（`readEGM2008` の既定出力と整合）。codegen 対応 (`%#codegen`)。

#### `readEGM2008`
```matlab
EGM = readEGM2008(EGM, deg, normalized)
```
EGM2008 係数ファイルを読み、構造体 `EGM` に係数行列を格納して返す。
- **入力**: `EGM` (struct) — 係数を書き込む構造体; `deg` — 読み込む次数; `normalized` — 省略/0 で非正規化係数、非0 で正規化係数（既定 0）
- **出力**: `EGM` (struct) — `.Cnm`, `.Snm`（各 (deg+1)×(deg+1)）を追加した構造体
- **Note**: 係数ファイル名は関数内で `'EGM2008_to2190_TideFree.txt'` にハードコードされており引数では渡さない。doc コメントおよび `egm2008.m` の使用例（`[EGM.Cnm, EGM.Snm] = readEGM2008('...txt', deg)`）は旧シグネチャで、現行の「構造体入力・構造体出力」と一致しない。

#### `shadow`
```matlab
nu = shadow(satI, sunI, rS, rE)
```
円錐影モデルによる地球影関数。
- **入力**: `satI` (nx3) — 慣性系の衛星位置; `sunI` (nx3) — 慣性系の太陽位置; `rS` — 太陽半径; `rE` — 地球半径
- **出力**: `nu` — 1: 日照、0: 本影（食）
- **Note**: 全変数の単位を統一すること（m でも km でも可）。

### 日心座標・定数

#### `earthVSOP87`
```matlab
[lon, lat, r] = earthVSOP87(jd, earthVSOP)
```
VSOP87 級数により、指定 Julian day における地球の日心黄経・黄緯・距離を計算する（当該日付の平均黄道・分点基準）。
- **入力**: `jd` (nx1, day) — Julian day; `earthVSOP` — `vsopConst` が返す係数テーブル
- **出力**: `lon` (mx1, rad) — 日心黄経; `lat` (mx1, rad) — 日心黄緯; `r` (mx1, AU) — 日心距離
- **Note**: doc は第2引数を `const` と記すが、実際は VSOP 係数 `earthVSOP` を渡す。

#### `vsopConst`
```matlab
earthVSOP = vsopConst
```
VSOP87 の係数テーブルを構造体 `earthVSOP` として返す（定数関数）。
- **出力**: `earthVSOP` — VSOP87 係数

#### `orbitConst`
```matlab
const = orbitConst()
```
軌道伝播に必要な物理定数を構造体 `const` にまとめて返す（定数関数）。
- **出力**: `const` — 地球・太陽・月の重力定数（`.GE/.GEm/.GEday/.GS/.GM/...`）、J2、地球半径 `.RE/.REm`、自転角速度 `.WE`、扁平率 `.fE`、J2000.0 平均黄道傾斜 `.EPS0`、AU、太陽定数 `.S0`、光速 `.c`、`.J2000` 等を含む

### 測地座標変換

#### `ecef2LatLonH`
```matlab
[lat, lon, h] = ecef2LatLonH(r, const)
```
ECEF 系の衛星位置を測地緯度・経度・高度へ変換する。
- **入力**: `r` (nx3, km) — ECEF 系位置; `const` (struct) — 軌道定数（地球半径・扁平率を使用）
- **出力**: `lat` (nx1, rad) — 測地緯度; `lon` (nx1, rad) — 経度; `h` — 楕円体高

#### `geodetic2Geocentric`
```matlab
r = geodetic2Geocentric(lat, lon, h, RE, f)
```
測地緯度・経度・高度から ECEF 位置ベクトルを計算する。
- **入力**: `lat` (nx1, rad) — 測地緯度; `lon` (nx1, rad) — 経度; `h` — 高度; `RE` — 地球赤道半径; `f` — 扁平率
- **出力**: `r` (nx3) — ECEF 系位置
- **Note**: `h`, `RE` の単位を揃えること。出力位置の単位は `RE`・`h` に一致する。

#### `geocentric2Geodetic`
```matlab
[lon, lat, h] = geocentric2Geodetic(x, y, z, a, f)
```
地心直交座標から測地緯度・経度・高度を反復計算する。
- **入力**: `x, y, z` — 地心直交座標; `a` — 参照楕円体の赤道半径; `f` — 扁平率 (a−b)/a
- **出力**: `lon` (rad, −π〜π) — 測地経度; `lat` (rad, −π/2〜π/2) — 測地緯度; `h` — 測地高（`a` と同単位）

### 恒星時・地球回転角

#### `gmst`
```matlab
GMST = gmst(jd)
```
Julian day から Greenwich 平均恒星時を計算する。
- **入力**: `jd` (nx1, day) — Julian day
- **出力**: `GMST` (nx1, rad) — Greenwich 平均恒星時

#### `gast`
```matlab
GAST = gast(jd, const)
```
IAU-76/FK5 に基づき Greenwich 視恒星時を計算する（章動による赤経差を含む）。
- **入力**: `jd` (nx1, day) — Julian day; `const` (struct) — 軌道定数
- **出力**: `GAST` (nx1, rad) — Greenwich 視恒星時

#### `era`
```matlab
theta = era(jdUT1)
```
IAU-2006/2000 (CIO ベース) の地球回転角を計算する。
- **入力**: `jdUT1` (nx1, day) — UT1 の Julian day
- **出力**: `theta` (rad) — 地球回転角

#### `earthW`
```matlab
w = earthW(lod)
```
LOD（1日の長さ）から地球自転角速度を計算する。
- **入力**: `lod` (nx1, s) — length of day
- **出力**: `w` (nx1, rad/s) — 自転角速度ノルム
- **Note**: 基準値 WE=7.292115146706979e-5 rad/s は `orbitConst().WE` と同期させること（コード内コメントで明記）。

#### `wobble`
```matlab
sPrime = wobble(tTT)
```
極運動の年周ゆらぎ s'（IAU-2006/2000, CIO ベース）を計算する。
- **入力**: `tTT` — TT のユリウス世紀
- **出力**: `sPrime` (rad) — 年周ゆらぎ角

### 歳差・章動

#### `obliquity`
```matlab
e = obliquity(jd)
```
指定元期における黄道傾斜角（平均）を計算する。
- **入力**: `jd` (nx1, day) — Julian day
- **出力**: `e` (rad) — 平均黄道傾斜角

#### `nutation`
```matlab
[dPsi, dEpsi] = nutation(jd, const)
```
平均赤道から真赤道 (TOD) への章動量（黄経・傾斜）を計算する。
- **入力**: `jd` (day) — Julian day; `const` (struct) — 軌道定数
- **出力**: `dPsi` (rad) — 黄経章動 dψ; `dEpsi` (rad) — 傾斜章動 dε

#### `nutationDCM`
```matlab
dcm = nutationDCM(jd, const)
```
MOD から TOD への章動回転行列を計算する (IAU-76/FK5)。
- **入力**: `jd` (day, scalar) — Julian day; `const` (struct) — 軌道定数
- **出力**: `dcm` (3x3) — MOD→TOD の DCM

#### `nutationQ`
```matlab
q = nutationQ(jd, scalar, const)
```
MOD から TOD への章動クォータニオンを計算する。
- **入力**: `jd` (nx1, day) — Julian day; `scalar` — クォータニオン定義（0: スカラー先頭 [q0,q1,q2,q3]、4: スカラー末尾 [q1,q2,q3,q4]）; `const` (struct) — 軌道定数
- **出力**: `q` (nx4) — MOD→TOD のクォータニオン

#### `precession`
```matlab
[zeta, z, theta, eta, Pi_, p] = precession(jd0, jd1, const)
```
平均から TOD への歳差角一式を計算する。
- **入力**: `jd0` (day, scalar) — 基準 Julian day; `jd1` (day, scalar) — 対象 Julian day; `const` (struct) — 軌道定数
- **出力**: `zeta, z, theta` (rad) — 歳差角; `eta` (rad) — 2黄道間の傾斜角; `Pi_` (rad) — 初期分点から黄道交点までの角; `p` (rad) — 黄経合成歳差

#### `precessionDCM`
```matlab
dcm = precessionDCM(jd0, jd1, const)
```
歳差の方向余弦行列を計算する（`jd0=J2000.0` のとき J2000(ICRF)→MOD の回転）。
- **入力**: `jd0` (day, scalar) — 基準 Julian day; `jd1` (day, scalar) — 対象 Julian day; `const` (struct) — 軌道定数
- **出力**: `dcm` (3x3) — jd0→jd1 の歳差 DCM

#### `precessionQ`
```matlab
q = precessionQ(jd0, jd1, scalar, const)
```
歳差クォータニオンを計算する（`jd0=J2000.0` のとき ICRF→TOD の回転）。
- **入力**: `jd0` (day, scalar) — 基準 Julian day; `jd1` (day, scalar) — 対象 Julian day; `scalar` — クォータニオン定義（0 または 4）; `const` (struct) — 軌道定数
- **出力**: `q` (nx4) — 歳差クォータニオン

#### `precessionNutation`
```matlab
dcm = precessionNutation(jdTT, iau06, dX, dY)
```
IAU-2006/2000 理論に基づき歳差＋章動の DCM（CIRS→GCRF）を計算する。
- **入力**: `jdTT` — TT のユリウス世紀; `iau06` — `readIAU06` が返す係数（X/Y/s テーブル）; `dX, dY` — 補正項（任意）
- **出力**: `dcm` (3x3) — CIRS→GCRF の DCM
- **Note**: `iau06x.dat, iau06y.dat, iau06z.dat` が必要（`readIAU06` で読み込む）。doc の入力欄は `iau06` を明記していないが第2引数として必須。

#### `earthFullRotQ`
```matlab
q = earthFullRotQ(jd0, jd1, scalar, const)
```
歳差と章動を含む地球回転クォータニオンを計算する（`nutationQ` と `precessionQ` を合成）。
- **入力**: `jd0` (day, scalar) — 基準 Julian day; `jd1` (nx1, day) — 対象 Julian day; `scalar` — クォータニオン定義（0 または 4）; `const` (struct) — 軌道定数
- **出力**: `q` (nx4) — jd0→jd1 の地球回転クォータニオン

#### `earthNutationPrecessionQ`
```matlab
q = earthNutationPrecessionQ(jd0, jd1, scalar, const)
```
歳差と章動を含む地球回転クォータニオンを計算する（極運動は含まない）。`earthFullRotQ` とタイトル・シグネチャがほぼ同一。
- **入力**: `jd0` (day, scalar) — 基準 Julian day; `jd1` (nx1, day) — 対象 Julian day; `scalar` — クォータニオン定義（0 または 4）; `const` (struct) — 軌道定数
- **出力**: `q` (nx4) — jd0→jd1 の地球回転クォータニオン
- **Note**: 極運動を含まない旨を note で明記。`earthFullRotQ` と実質同機能で、どちらを使うか要確認。

### フレーム変換 (ITRF/GCRF/PEF)

#### `itrf2gcrf`
```matlab
dcm = itrf2gcrf(jd, EOP)
```
IAU-2006/2000 (CIO アプローチ) による ITRF→GCRF 回転行列を計算する（位置座標のみ）。
- **入力**: `jd` (day) — Julian day; `EOP` (struct) — `readEOP` の出力（`.dataAll`, `.iau06`, `.leapJD` を含む）
- **出力**: `dcm` (3x3) — ITRF→GCRF の回転行列
- **Note**: 位置座標の変換のみ。速度変換には地球自転速度の考慮が別途必要。codegen 対応 (`%#codegen`)。

#### `qITRF2gcrf`
```matlab
q = qITRF2gcrf(scalar, jd)
```
`itrf2gcrf` のラッパで、ITRF→GCRF クォータニオンを引数順を変えて返す。
- **入力**: `scalar` — クォータニオン定義（0 または 4）; `jd` (day) — Julian day
- **出力**: `q` (1x4) — ITRF→GCRF クォータニオン
- **Note**: 位置座標の変換のみ。

#### `pef2itrf`
```matlab
R = pef2itrf(xp, yp)
```
極運動パラメータから PEF→ITRF の DCM を計算する。
- **入力**: `xp, yp` (rad) — 地球姿勢パラメータ（極運動）
- **出力**: `R` (3x3) — PEF→ITRF の DCM
- **Note**: Vallado の教科書は ITRF→PEF だが、本関数はその逆（PEF→ITRF）である点に注意。

### フレーム変換 (軌道要素ベース)

#### `mod2J2000`
```matlab
[iJ, OmeJ, wJ] = mod2J2000(jd, i, Ome, w, const)
```
MOD 系の軌道要素（i, Ω, w）を J2000 系へ変換する（Newcomb の歳差定数）。
- **入力**: `jd` (day) — Julian day; `i` (rad) — 軌道傾斜角; `Ome` (rad) — 昇交点経度; `w` (rad) — 近点引数; `const` (struct) — 軌道定数
- **出力**: `iJ` (rad) — J2000.0 の傾斜角; `OmeJ` (rad) — J2000.0 の昇交点経度; `wJ` (rad) — J2000.0 の近点引数

#### `teme2J2000`
```matlab
[iJ, OmeJ, wJ] = teme2J2000(jd, i, Ome, w, const)
```
TEME 系の軌道要素を J2000(FK5) 系へ変換する（Newcomb の歳差定数）。
- **入力**: `jd` (day) — Julian day; `i` (rad) — 軌道傾斜角; `Ome` (rad) — 昇交点経度; `w` (rad) — 近点引数; `const` (struct) — 軌道定数
- **出力**: `iJ` (rad) — J2000.0 の傾斜角; `OmeJ` (rad) — J2000.0 の昇交点経度; `wJ` (rad) — J2000.0 の近点引数
- **Note**: doc の入力欄には `const` の記載がないが、シグネチャは第5引数として `const` を要求する。

#### `tod2Mod`
```matlab
[iMod, OmeMod, wMod] = tod2Mod(jd, i, Ome, w, const)
```
TOD 系の軌道要素を MOD 系へ変換する（IAU-1980 章動理論）。
- **入力**: `jd` (day) — Julian day; `i` (rad) — 軌道傾斜角; `Ome` (rad) — 昇交点経度; `w` (rad) — 近点引数; `const` (struct) — 軌道定数
- **出力**: `iMod` (rad) — MOD の傾斜角; `OmeMod` (rad) — MOD の昇交点経度; `wMod` (rad) — MOD の近点引数

### データ読み込み (EOP/IAU06/JB2008/TLE)

#### `eop`
```matlab
output = eop(year, month, day, eopDataAll)
```
指定日付の地球姿勢パラメータ (EOP) を EOP データ全体から抽出する。
- **入力**: `year, month, day` (各 nx1) — 対象日付; `eopDataAll` (struct) — 読み込み済み EOP データ全体
- **出力**: `output` (struct) — `.mjd`（修正 Julian day）、`.xp/.yp`（極運動, rad）、`.dUT1`（ΔUT1）、`.dX/.dY`（rad）、`.lod`（length of day）
- **Note**: IERS 由来のコメント除去済みデータが必要（サンプル: `EOP_20_C04_one_file_1962-now.txt`）。

#### `readEOP`
```matlab
EOP = readEOP(fName)
```
EOP データファイルを読み込む。
- **入力**: `fName` — EOP データファイル名
- **出力**: `EOP` (struct) — EOP データ

#### `readIAU06`
```matlab
iau06 = readIAU06(~)
```
IAU-2006/2000 の歳差・章動計算用 IAU06 (x,y,z) 係数を読み込む。
- **入力**: なし（プレースホルダ引数）
- **出力**: `iau06` — IAU06 係数
- **Note**: `iau06x.dat, iau06y.dat, iau06z.dat` が必要。

#### `readJB2008`
```matlab
[PC, EOP, SOL, DTC] = readJB2008
```
Jacchia–Bowman 2008 大気モデルの係数群を読み込む。
- **入力**: なし
- **出力**: `PC` — 係数; `EOP` — 地球姿勢パラメータ; `SOL` — 宇宙天気データ; `DTC` — 地磁気嵐 DTC 値

#### `readTLE`
```matlab
tle = readTLE(tleName, const, tool)
```
TLE（2 行軌道要素）ファイルを読み、平均軌道要素・元期・識別子等を構造体で返す。
- **入力**: `tleName` — TLE ファイル名; `const` (struct) — 軌道定数; `tool` — 時刻変換に使うツール `'yoshimuLibrary'`（既定）/`'SPICE'`/`'MATLAB'`
- **出力**: `tle` (struct) — `.satName/.satID/.launchYear/.launchNum/.launchPiece`、`.oe`（平均軌道要素 a[km], e, i, Ω, w, M[rad]）、`.jd`（TT の Julian day）、`.n`（平均運動 rev/day）等
- **Note**: doc は SPICE 利用フラグを `spiceFlag`（on/off）と記すが、現行シグネチャは第3引数 `tool` で `{'yoshimuLibrary','SPICE','MATLAB'}` を選ぶ列挙型。

### Scripts

#### `compareLegendre`
```matlab
compareLegendre   % スクリプト
```
随伴 Legendre 関数の 2 実装（MATLAB 組み込み `legendre` を使う方法と再帰計算）を計時比較するスクリプト。次数が非常に大きい場合（〜150 以上）を除き組み込み関数側が速い、というメモが冒頭にある。関数ではなくトップレベルのスクリプト。

---

## 初期軌道決定 (`orbitDetermination/`)

3 時点の観測（角度観測または位置ベクトル）から軌道を初期決定する古典的アルゴリズム（Gauss 法・Double-r 法・Gibbs 法）を収める。いずれも観測は時系列順 `t1 < t2 < t3` に並んでいることを前提とし、位置・速度は行方向にサンプル（時点）、列方向に成分 (x, y, z) を並べる横ベクトル規約で扱う。角度観測系（`gauss` / `doubleR`）は方位・仰角を `aziele` (3x2)、観測点の慣性系位置を (3x3) で受け取り、位置ベクトル系（`gibbs`）は 3 本の 1x3 位置ベクトルを直接受け取る。重力定数 `mu` は入力位置ベクトルと単位系を揃えること（km を使うなら km^3/s^2）。実行例は本ライブラリの `examples/orbitDetermination/`（`exampleGauss` / `exampleGibbs` / `exampleDoubleR`）にある。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `gauss` | Gauss 法による角度観測 3 点からの初期軌道決定（8 次方程式の初期値を figure 上のクリックで対話的に選ぶ） |
| `doubleR` | Double-r 法による角度観測 3 点からの初期軌道決定（2 変数 Newton 反復） |
| `gibbs` | Gibbs 法による位置ベクトル 3 点からの速度決定 |

#### `gauss`
```matlab
[r2, v2] = gauss(t, aziele, obsECI, rRange, mu)
```
Gauss 法。3 時点の角度観測（方位・仰角）と観測点位置から中央時点の位置・速度ベクトルを求める。8 次方程式の求根は、`rRange` 上で評価した多項式を figure に描き、ユーザがライン上をクリックして初期値を与え、そこから Newton 法で収束させる（対話的処理を含む）。
- **入力**: `t` (3x1, 日 [day]) — 観測時刻; `aziele` (3x2, rad) — 各観測の方位角・仰角（トポセントリック赤道座標系）; `obsECI` (3x3) — 各観測時の観測点位置ベクトル（慣性系、行=時点）; `rRange` (1xN) — 8 次方程式の求根範囲兼プロット横軸（初期値の候補格子）; `mu` (1x1) — 地球重力定数
- **出力**: `r2` (1x3) — 中央時点 t2 の位置ベクトル（慣性系）; `v2` (1x3) — 中央時点 t2 の速度ベクトル
- **Note**: `waitforbuttonpress` と figure クリックで 8 次方程式の初期解を選ぶため、バッチ／非対話実行には不向き。`r2` は出力（位置ベクトル）であり、内部で観測点 `R2` などとは別物である点に注意。
- *See also*: `gibbs`

#### `doubleR`
```matlab
[r2, v2] = doubleR(t, aziele, rObs, mu, RE)
```
Double-r 法。2 つの動径距離 r1, r2 を未知数とし、2 変数の Newton 反復（差分による偏微分）で観測を満たす軌道を求める。中央時点の位置・速度ベクトルを返す。
- **入力**: `t` (3x1, 日 [day]) — 観測時刻（内部で `day2s` により秒へ換算）; `aziele` (3x2, rad) — 各観測の方位角・仰角; `rObs` (3x3) — 観測点位置ベクトル（慣性系、行=時点）; `mu` (1x1) — 地球重力定数; `RE` (1x1) — 地球半径。反復の収束判定 (`TOL = 1e-8 * RE`) と初期推定 (`r1 = 2.0*RE`, `r2 = 2.01*RE`) のスケールに使う
- **出力**: `r2` (1x3) — 中央時点 t2 の位置ベクトル（慣性系）; `v2` (1x3) — 中央時点 t2 の速度ベクトル
- **Note**: ヘッダコメントは入力に `rRange` を挙げているが、実際のシグネチャに `rRange` は無く、代わりに未文書の `RE`（地球半径）を取る。単位系は `mu`・`RE` とベクトルで揃えること。内部で呼ぶ `day2s`（日→秒換算）はこのライブラリに存在せず、別途パス上に必要（`conversion/` にあるのは逆変換の `s2day` のみ）。
- *See also*: `gibbs`

#### `gibbs`
```matlab
[v1, v2, v3] = gibbs(r1, r2, r3, mu)
```
Gibbs 法。同一軌道上の 3 つの位置ベクトルから、各時点の速度ベクトルを解析的に求める。3 点が同一平面上にない場合はエラーを送出する。
- **入力**: `r1, r2, r3` (各 1x3) — t1 < t2 < t3 の順に並んだ地心位置ベクトル; `mu` (1x1) — 重力定数（位置ベクトルと単位を揃える。位置が km なら km^3/s^2）
- **出力**: `v1, v2, v3` (各 1x3) — 各時点の速度ベクトル
- **Note**: 位置ベクトルは時系列順 (t1 < t2 < t3) であること。共面性チェックは `u1 · (r2×r3) > 1e-3` で非共面と判定しエラーとなる（単位ベクトル u1 との内積に対する固定しきい値）。
- *See also*: `orbitConst`

---

## 相対軌道 (ROE / 視線角) (`relativeOrbit/`)

相対軌道要素 (Relative Orbital Elements: ROE) と相対運動を扱う関数群。chief/deputy の絶対軌道要素 (nx6: a, e, i, Ω, w, f または M) ↔ ROE (nx6: δa, δλ, δex, δey, δix, δiy) の相互変換、ROE から RTN (Radial-Transverse-Normal) 座標系での相対位置・速度へのマッピング、および角度のみ航法 (angles-only navigation) 用の視線角 (LOS: azimuth/elevation) の計算を提供する。行方向=時系列サンプル、列方向=状態要素の横ベクトル規約 (nx6, nx3, nx1) に従う。多くの関数が `anomalyFlag`/`flag` (1=真近点角, 0=平均近点角) と重力定数 `GE` を取り、`GE` と位置・速度・semi-major axis の単位 (km か m) を揃える必要がある。ROE 定義は δα=(a_d-a)/a, δλ=(u_d-u)+(Ω_d-Ω)cos i, δex=e_xd-e_x, δey=e_yd-e_y, δix=i_d-i, δiy=(Ω_d-Ω)sin i。`calcRelPosVelAtti` のみ `chief`/`deputy` 構造体と `const` 構造体を受け取る。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `oe2roe` | chief/deputy の絶対軌道要素から相対軌道要素 (ROE) を計算 |
| `roe2DeputyOE` | ROE と chief の軌道要素から deputy の絶対軌道要素を復元 |
| `roe2rtn` | ROE を RTN 座標系での相対位置・速度に (1次) マッピング |
| `oe2los` | 絶対軌道要素から厳密な相対位置を求め視線角 (azimuth/elevation) を計算 |
| `roe2mappedLOS` | ROE の1次 RTN マッピングを介した近似視線角を計算 |
| `calcRelPosVelAtti` | chief/deputy 構造体から相対位置・速度・姿勢を一括計算 |

### Core (ROE 変換)

#### `oe2roe`
```matlab
roe = oe2roe(chief, deputy, anomalyFlag)
```
chief と deputy の絶対軌道要素から相対軌道要素 (ROE) を計算する。
- **入力**: `chief` (nx6, [km, -, rad, rad, rad, rad]) — chief の絶対軌道要素 [a, e, i, Ω, w, f または M]; `deputy` (nx6, 同上) — deputy の絶対軌道要素; `anomalyFlag` (1x1) — 1=真近点角, 0=平均近点角
- **出力**: `roe` (nx6, [-, rad, -, -, rad, rad]) — [δa, δλ, δex, δey, δix, δiy]
- **Note**: `δλ` は wrapPi で [-π, π) に折り返される。RAAN/w/f(M) は内部で mod 2π される。
- *See also*: `roe2DeputyOE`

#### `roe2DeputyOE`
```matlab
deputyOE = roe2DeputyOE(roe, chiefOE, anomalyFlag)
```
ROE と chief の絶対軌道要素から deputy の絶対軌道要素を復元する (`oe2roe` の逆変換)。
- **入力**: `roe` (nx6) — 相対軌道要素 [δa, δλ, δex, δey, δix, δiy]; `chiefOE` (nx6, [km, -, rad, rad, rad, rad]) — chief の絶対軌道要素; `anomalyFlag` (1x1) — 1=真近点角, 0=平均近点角
- **出力**: `deputyOE` (nx6, [km, -, rad, rad, rad, rad]) — deputy の絶対軌道要素。`anomalyFlag`=1 なら6列目は真近点角 f、0 なら平均近点角 M
- **Note**: chief の傾斜角 iC が eps 未満 (near-equatorial) のときは δiy/sin(iC) の特異点を避けるため RAAN を chief と同一にする。出典 Di Mauro (2019, JGCD)。
- *See also*: `oe2roe`

#### `roe2rtn`
```matlab
[xRTN, vRTN] = roe2rtn(roe, chiefOE, flag, GE)
```
ROE を chief 中心の RTN 座標系での相対位置・速度に1次 (線形) マッピングする。
- **入力**: `roe` (nx6) — 相対軌道要素; `chiefOE` (nx6, [km または m, -, rad, rad, rad, rad]) — chief の絶対軌道要素; `flag` (1x1) — 1=真近点角, 0=平均近点角; `GE` (1x1) — 地球重力定数 (単位は位置・速度と統一)
- **出力**: `xRTN` (nx3, km または m) — deputy の RTN 相対位置; `vRTN` (nx3, km/s または m/s) — deputy の RTN 相対速度
- **Note**: `GE` と semi-major axis (chiefOE 1列目) の単位を揃えること。円軌道近似の1次マッピング。

### Variants (視線角 / LOS)

#### `oe2los`
```matlab
[azi, ele] = oe2los(chiefOE, deputyOE, anomalyFlag, GE, chiefQbi)
```
chief/deputy の絶対軌道要素から厳密な相対位置を求め、角度のみ航法用の視線角を計算する。`chiefQbi` 省略時は RTN 座標系、指定時は chief の機体固定系で LOS を評価する。
- **入力**: `chiefOE` (nx6, [km または m, -, rad, rad, rad, rad]) — chief の絶対軌道要素; `deputyOE` (nx6, 同上) — deputy の絶対軌道要素; `anomalyFlag` (1x1) — 1=真近点角, 0=平均近点角; `GE` (1x1) — 地球重力定数 (位置・速度と単位統一); `chiefQbi` (省略可, クォータニオン) — chief の慣性系→機体系クォータニオン。省略時は RTN 系で LOS を計算
- **出力**: `azi` (nx1, rad) — R-T 平面内の R 軸からの方位角 (atan2(rel_T, rel_R)); `ele` (nx1, rad) — R-T 平面から N 方向への仰角
- **Note**: `GE` と位置・速度の単位を揃えること。出典 Sullivan et al. (Generalized Angles-Only Navigation, JGCD)。`arguments` ブロックは `chiefQbi (1,4)` と宣言するが、`qRotation` で nx3 の相対ベクトルに適用されるため時系列入力では nx4 が想定される点に注意。
- *See also*: (doc は `roe2losApprox` を挙げるが未存在。近似版は `roe2mappedLOS`)

#### `roe2mappedLOS`
```matlab
[azi, ele] = roe2mappedLOS(roe, chiefOE, flag, GE, chiefQbo)
```
ROE を1次 RTN マッピング (`roe2rtn`) で相対位置に変換してから近似視線角を計算する (`oe2los` の近似版)。`chiefQbo` 省略時は RTN 系、指定時は chief の機体固定系で評価。
- **入力**: `roe` (nx6) — 相対軌道要素; `chiefOE` (nx6, [km または m, -, rad, rad, rad, rad]) — chief の絶対軌道要素; `flag` (1x1) — 1=真近点角, 0=平均近点角; `GE` (1x1) — 地球重力定数 (semi-major axis と単位統一); `chiefQbo` (省略可, クォータニオン) — chief の RTN 系→機体系クォータニオン。省略時は RTN 系
- **出力**: `azi` (nx1, rad) — R-T 平面内の R 軸からの方位角 (atan2(rel_T, rel_R)); `ele` (nx1, rad) — R-T 平面から N 方向への仰角
- **Note**: `GE` と semi-major axis の単位を揃えること。1次 ROE マッピングを介するため厳密な `oe2los` に対し近似。出典 Di Mauro (2019, JGCD)。
- *See also*: `oe2roe`

### Helpers (構造体一括計算)

#### `calcRelPosVelAtti`
```matlab
[chief, deputy, rel] = calcRelPosVelAtti(chief, deputy, anomalyFlag, const)
```
chief/deputy 構造体から慣性系の絶対位置・速度、RTN 系での相対位置・速度 (非線形/ROE 両方)、および (姿勢が与えられていれば) 相対姿勢・角速度を一括計算する。
- **入力**: `chief` (構造体) — フィールド `oe` (nx6 軌道要素) 必須、任意で `q` (慣性系姿勢クォータニオン), `n` (平均運動), `w` (角速度); `deputy` (構造体) — 同上のフィールド; `anomalyFlag` (1x1) — 1=真近点角, 0=平均近点角; `const` (構造体) — フィールド `GE` (地球重力定数) を使用
- **出力**: `chief` — `rI`/`vI` (慣性系位置・速度), `q` あるとき `qoi`/`qbo` を追加した構造体; `deputy` — `rI`/`vI`, `q` あるとき `qbo`/`woi`/`wbo` を追加した構造体; `rel` — 相対量構造体 (`rNonlinI`/`rNonlinRTN`/`vNonlinI`/`vNonlinRTN` の非線形相対位置速度, `roe`, `rMappedRTN` の ROE マッピング相対位置)
- **Note**: 姿勢・角速度の計算は `chief` に `q` フィールドがある場合のみ実行される。`oe2rv`, `oe2roe`, `roe2rtn`, `triad` 等の外部関数に依存。ROE の δλ は wrapPi で [-π, π) に折り返される。

---

## 太陽・月の位置と重力 (`sunMoon/`)

この節では、太陽と月の位置（J2000.0 慣性系）および人工衛星に作用する第三体重力（摂動加速度）を計算する関数群を扱う。位置は中精度モデル（太陽 = VSOP87、月 = ELP2000-82 理論）で求め、必要なら SPICE/MICE カーネルによる高精度計算に切り替えられる。時刻・状態は行方向を時系列サンプル（n サンプル）とする横ベクトル規約（位置は n×3）で受け渡し、`const` 構造体（`orbitConst` が生成、`EPS0` / `GM` / `GS` / `J2000` 等を保持）と VSOP/ELP 係数構造体（`earthVSOP87` / `vsopConst`、`readELP`）に依存する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `sun` | 太陽の J2000.0 慣性系位置ベクトルを計算 |
| `sunLonLatR` | その日付の平均黄道・分点基準での太陽の幾何黄経・黄緯・距離 |
| `sunG` | 太陽による第三体重力加速度（衛星に働く摂動） |
| `moon` | 月の J2000.0 慣性系位置ベクトルを計算 |
| `moonLonLatR` | J2000.0 の平均黄道・分点基準での月の幾何黄経・黄緯・距離 |
| `moonELP` | ELP2000-82 理論による月の幾何黄経・黄緯・距離（date 基準） |
| `moonG` | 月による第三体重力加速度（衛星に働く摂動） |
| `readELP` | 月の ELP 係数を CSV から読み込み `ELP` 構造体を生成 |
| `makeELPcoeff` | ELP 係数表を CSV（`ELPcoeffA.csv` / `ELPcoeffB.csv`）に書き出す (スクリプト) |

### 太陽 (Sun)

#### `sun`
```matlab
sunPos = sun(jd, const, earthVSOP)
```
指定ユリウス日における太陽位置を J2000.0 慣性系で返す。`sunLonLatR` で得た黄経・黄緯・距離から直交座標を組み立て、平均黄道傾斜角 `-EPS0` の X 軸回転（クォータニオン）で J2000.0 系に変換する。
- **入力**: `jd` (:,1, day) — UTC のユリウス日（列ベクトル化される）; `const` — 軌道定数構造体（`EPS0` を使用）; `earthVSOP` — 地球 VSOP87 係数
- **出力**: `sunPos` (n×3, km) — J2000.0 慣性系での太陽位置ベクトル
- *See also*: `sunLonLatR`, `vsopConst`

#### `sunLonLatR`
```matlab
[lon, lat, r] = sunLonLatR(jd, const, earthVSOP)
```
その日付の平均黄道・分点を基準とした太陽の幾何黄経・黄緯・距離を計算する。`earthVSOP87` で得た地球の日心座標に黄経 +180°・黄緯符号反転を施して太陽の地心座標に変換し、歳差補正を掛ける。
- **入力**: `jd` (:,1, day) — ユリウス日; `const` — 軌道定数構造体; `earthVSOP` — 地球 VSOP87 係数
- **出力**: `lon` (n, rad) — 太陽の地心黄経; `lat` (n, rad) — 太陽の地心黄緯; `r` (n, AU) — 太陽の地心距離
- **Note**: 距離 `r` の単位は AU（他の位置関数の km とは異なる）。`sun` 側で `au2km` により km へ換算される。
- *See also*: `orbitConst`, `vsopConst`, `precession`

#### `sunG`
```matlab
[aSun, sunI] = sunG(jd, rVec, const, earthVSOP, options)
```
太陽による第三体重力（直接項 − 間接項）の加速度を慣性系で計算する。既定では `sun` で位置を求め、`options.SPICE='on'` の場合は SPICE（`cspice_spkpos`）で太陽位置を取得する。
- **入力**: `jd` (:,1, day) — ユリウス日; `rVec` (n×3, km) — 慣性系での衛星位置; `const` — 軌道定数構造体（`GS` を使用）; `earthVSOP` — 地球 VSOP87 係数; `options.SPICE` (char) — `'on'`/`'off'`（既定 `'off'`）
- **出力**: `aSun` (n×3, km/s^2) — 太陽の重力加速度; `sunI` (n×3, km) — 慣性系（GCRF）での太陽位置
- **Note**: `SPICE='on'` には SPICE/MICE カーネルのロードが必要。
- *See also*: `orbitConst`, `vsopConst`, `precession`

### 月 (Moon)

#### `moon`
```matlab
moonPos = moon(jd, const, ELP)
```
指定ユリウス日における月位置を J2000.0 慣性系で返す。`moonLonLatR` で得た黄経・黄緯・距離から直交座標を組み立て、`-EPS0` の X 軸回転（クォータニオン）で J2000.0 系に変換する。
- **入力**: `jd` (:,1, day) — ユリウス日; `const` — 軌道定数構造体（`EPS0` を使用）; `ELP` — ELP 係数構造体（`readELP` で生成）
- **出力**: `moonPos` (n×3, km) — J2000.0 慣性系での月位置ベクトル
- *See also*: `sun`

#### `moonLonLatR`
```matlab
[lonM, latM, rM] = moonLonLatR(jd, const, ELP)
```
J2000.0 の平均黄道・分点を基準とした月の幾何黄経・黄緯・距離を計算する。`moonELP`（利用可能なら codegen 済み `moonELP_mex`）で date 基準の値を求め、歳差補正で J2000.0 基準に変換する。
- **入力**: `jd` (:,1, day) — ユリウス日; `const` — 軌道定数構造体（`J2000` を使用）; `ELP` — ELP 係数構造体
- **出力**: `lonM` (n, rad) — 月の地心黄経; `latM` (n, rad) — 月の地心黄緯; `rM` (n, km) — 月の地心距離
- *See also*: `sun`

#### `moonELP`
```matlab
[lon, lat, r] = moonELP(jd, ELP)
```
ELP2000-82 理論（Meeus の簡略級数）で月の幾何黄経・黄緯・距離を計算する。その日付（date）の平均黄道・分点が基準。係数表の各項を D/M/M'/F の引数で級数展開し、太陽離心率補正 E を掛けて総和する。
- **入力**: `jd` (:,1, day) — ユリウス日; `ELP` — ELP 係数構造体（`ELP.a` = 黄経・距離表、`ELP.b` = 黄緯表）
- **出力**: `lon` (n, rad) — 月の地心黄経; `lat` (n, rad) — 月の地心黄緯; `r` (n, km) — 月の地心距離
- **Note**: 出力は date 基準。J2000.0 基準が必要なら `moonLonLatR` を使う。`%#codegen` 対応関数。
- *See also*: `sun`

#### `moonG`
```matlab
[aMoon, moonIJK] = moonG(jd, rVec, const, ELP, options)
```
月による第三体重力（直接項 − 間接項）の加速度を慣性系で計算する。既定では `moon` で位置を求め、`options.SPICE='on'` の場合は SPICE（`cspice_spkpos`）で月位置を取得する。
- **入力**: `jd` (:,1, day) — ユリウス日; `rVec` (n×3, km) — 慣性系での衛星位置; `const` — 軌道定数構造体（`GM` を使用）; `ELP` — ELP 係数構造体（SPICE 使用時は不要）; `options.SPICE` (char) — `'on'`/`'off'`（既定 `'off'`）
- **出力**: `aMoon` (n×3, km/s^2) — 月の重力加速度; `moonIJK` (n×3, km) — 慣性系での月位置
- **Note**: `ELP` は SPICE 使用時には不要。`SPICE='on'` には SPICE/MICE カーネルのロードが必要。
- *See also*: `orbitConst`, `sunG`

### 係数の入出力 (Helpers / Scripts)

#### `readELP`
```matlab
ELP = readELP(fNames)
```
月の ELP 係数を CSV から読み込み、`moon` / `moonELP` が使う `ELP` 構造体（`ELP.a`, `ELP.b`）を生成する。
- **入力**: `fNames` (1×2 cell, optional) — 係数 CSV のファイル名 2 つ。省略時は `{'ELPcoeffA.csv','ELPcoeffB.csv'}` を既定使用
- **出力**: `ELP` — `ELP.a`（黄経・距離表）と `ELP.b`（黄緯表）を持つ構造体
- *See also*: `moonG`

#### `makeELPcoeff` (スクリプト)
```matlab
makeELPcoeff   % スクリプト（引数・戻り値なし）
```
Meeus の Table 47.A / 47.B を埋め込んだ配列から、月位置計算用の ELP 係数 CSV（`ELPcoeffA.csv`, `ELPcoeffB.csv`）をカレントフォルダに書き出す。`readELP` が読む係数ファイルの生成用。
- **Note**: 冒頭で `clc/clear/cls` を実行し、`writematrix` でカレントフォルダに CSV を出力する（既存ファイルを上書き）。単発の準備スクリプトで、通常の計算では `readELP` 経由で使う。

---

## 宇宙環境モデル (IGRF・大気密度) (`environment/`)

このディレクトリは、地球磁場 (IGRF-12) と超高層大気密度 (Jacchia–Roberts 1971 / Jacchia–Bowman 2008) を評価する環境モデルと、それらが必要とする太陽・地磁気指数 (CelesTrak SW-All.csv) の読み込みユーティリティをまとめる。IGRF 関数は緯度・経度をラジアンで受け取り (内部で deg 変換)、磁場を NED 座標 (北向き・東向き・下向き) の 1x3 横ベクトル [nT] で返す。大気密度モデルは高度の単位が関数ごとに異なる (jr1971 は [m]、jacciaBowman は [km]) 点に注意。太陽指数は `orbitConst` 由来の `const` 構造体や CelesTrak の `sw` 構造体に依存する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `igrf12` | IGRF-12 モデルで磁場ベクトルを計算 (mex 使用、単一時刻) |
| `geodeticIGRF` | 測地座標入力で IGRF 磁場を計算するラッパ (スカラ入力用) |
| `jr1971` | Jacchia–Roberts 1971 大気モデルで密度・温度・各種数密度を計算 |
| `jacciaBowman` | Jacchia–Bowman 2008 大気モデルで密度・温度を計算するラッパ |
| `loadSpaceWeather` | CelesTrak SW-All.csv を読み込み Kp/Ap/F10.7 構造体を返す |
| `lookupSolarGeoIndex` | 指定 JD の F10.7・F10.7a・Kp を取得 (JR1971 用) |
| `disturbances` | 地球・J2・太陽・月による摂動加速度を高度に対してプロット (スクリプト) |

### 地磁場モデル (IGRF)

#### `igrf12`
```matlab
b = igrf12(year, alt, lat, lon)
```
IGRF-12 モデルで指定時刻・地点の地球磁場ベクトルを計算する。`igrfsyn12.c` から生成した mex ファイルを呼ぶ。
- **入力**: `year` (1x1, 年) — 小数年 (例 2019.5); `alt` (1x1, km) — 海抜高度; `lat` (1x1, rad) — 北緯; `lon` (1x1, rad) — 東経
- **出力**: `b` (1x3, nT) — NED 座標系の磁場ベクトル [北向き, 東向き, 下向き]
- **Note**: 実行には `igrfsyn12` の mex ファイルが path 上に必要。`environment/` にはソース (`igrfsyn12.c`/`.h`) のみが置かれ、コンパイル済み mex は `examples/attitudeOrbit/igrfsyn12.mexmaci64` にある。緯度・経度は関数内部で deg に変換される。

#### `geodeticIGRF`
```matlab
b = geodeticIGRF(jd, lat, lon, alt, coefs)
```
Julian day と測地座標 (緯度・経度・高度) を入力に IGRF 磁場を計算するラッパ。スカラ入力を想定。
- **入力**: `jd` (1x1, day) — Julian day; `lat` (1x1, rad) — 測地緯度; `lon` (1x1, rad) — 測地経度; `alt` (1x1, km) — 測地高度; `coefs` — IGRF 係数構造体 (`igrfcoefs.mat` 相当)
- **出力**: `b` (1x3, nT) — NED 座標系の磁場ベクトル [北向き, 東向き, 下向き]
- **Note**: 内部で `jd2gc` により JD を暦日に変換し `igrfs(...)` を呼ぶが、**`igrfs` はライブラリ内に存在しない** (下記 不整合を参照)。緯度・経度は内部で deg に変換される。

### 大気密度モデル

#### `jr1971`
```matlab
result = jr1971(jd, phi_gd, lambda, h, F10, F10a, Kp)
```
Jacchia–Roberts 1971 大気モデルで、指定時刻・地点・高度における全大気密度・温度・主要 6 化学種の数密度を計算する。密度計算のロジックはこのファイル内のローカル関数群で完結している。
- **入力**: `jd` (1x1, day) — Julian day; `phi_gd` (1x1, rad) — 測地緯度; `lambda` (1x1, rad) — 経度; `h` (1x1, m) — 高度 (**メートル**で渡す); `F10` (1x1, sfu) — 10.7cm 太陽フラックス; `F10a` (1x1, sfu) — 81 日中央平均 F10.7; `Kp` (1x1) — Kp 地磁気指数 (3 時間遅延)
- **出力**: `result` (struct) — フィールド: `total_density` [kg/m^3], `temperature` [K], `exospheric_temperature` [K], `N2_number_density` / `O2_number_density` / `O_number_density` / `Ar_number_density` / `He_number_density` / `H_number_density` [1/m^3]
- **Note**: 高度 `h` は [m] で渡す (関数冒頭で km に変換)。`F10`/`F10a`/`Kp` は `lookupSolarGeoIndex(jd, 'SW-All.csv')` から取得できる。定数は `persistent` でキャッシュされる。

#### `jacciaBowman`
```matlab
[temp, rho] = jacciaBowman(jd, lon, lat, h, const, JB)
```
Jacchia–Bowman 2008 (JB2008) 大気モデルで温度と大気密度を計算するラッパ。太陽指数の遅延 (F10/S10 は 1 日、M10 は 2 日、Y10 は 5 日) や地磁気補正・座標変換を組み立てて `JB2008(...)` を呼ぶ。
- **入力**: `jd` (1x1, day) — Julian day; `lon` (1x1, rad) — 経度 (地心経度); `lat` (1x1, rad) — 地心緯度; `h` (1x1, km) — 測地高度; `const` — 軌道定数構造体 (`orbitConst` 由来); `JB` — JB2008 係数・EOP/太陽指数データ構造体 (`SOLdata`/`DTCdata`/`EOPdata` 等)
- **出力**: `temp` (1x1, K) — 温度; `rho` (1x1, kg/m^3) — 大気密度
- **Note**: 高度 `h` は [km] で渡す (jr1971 の [m] と異なる)。核となる `JB2008`・`finddays`・`IERS`・`JPL_Eph_DE430` 等はこのライブラリに含まれず、外部ツールボックス (Mahooti HPOP 系) が path 上に必要。第 2 引数はコード上「経度 (longitude)」として使われるが、doc コメントでは `alp: right ascension` と記載されており不一致 (下記)。
- *See also*: `jd2mjd`, `jd2gc`

### 太陽・地磁気指数ユーティリティ

#### `loadSpaceWeather`
```matlab
sw = loadSpaceWeather(csvPath)
```
CelesTrak の宇宙天気ファイル SW-All.csv を読み込み、Kp/Ap/F10.7 を含む構造体を返す。日付スパンと行数の連続性をチェックする。
- **入力**: `csvPath` — SW-All.csv のファイルパス
- **出力**: `sw` (struct) — `jdNoonStart`/`jdNoonEnd` (先頭・末尾行の正午 JD, 整数), `nDays` (行数=連続日数), `Kp` (Nx8, 3 時間ブロック実数), `Ap` (Nx8), `ApAvg` (Nx1), `F107obs` (Nx1) [sfu], `F107obsCenter81` (Nx1, 81 日中央平均観測値) [sfu]
- **Note**: CelesTrak の Kp は 10 倍で格納されているため、読み込み時に 10 で割って実数化する。データ取得元: https://celestrak.org/SpaceData/SW-All.csv

#### `lookupSolarGeoIndex`
```matlab
[F10, F10a, Kp] = lookupSolarGeoIndex(jd, sw)
```
指定 Julian day に対する太陽フラックス F10.7 (1 日ラグ)・81 日中央平均・Kp (ラグなし) を取得する。`jr1971` への入力を用意する用途。
- **入力**: `jd` (1x1, day) — Julian day (UTC); `sw` — `loadSpaceWeather` が返す構造体 (推奨) または SW-All.csv のパス文字列 (パス指定時は `persistent` でキャッシュ)
- **出力**: `F10` (1x1, sfu) — 10.7cm 太陽フラックス (1 日ラグ, 観測値); `F10a` (1x1, sfu) — 81 日中央平均 F10.7 (観測値); `Kp` (1x1) — jd 時刻の 3 時間ブロックの Kp (ラグなし)
- **Note**: F10/F10a は「前日正午」の値 (`floor(jd-0.5)`)、Kp は「当日」の 3 時間ブロック値 (`floor(jd+0.5)` と 時刻から `floor(hr/3)+1`) を参照する。範囲外・NaN の場合はエラーを送出する。

### スクリプト

#### `disturbances`
```matlab
disturbances   % スクリプト (引数なし)
```
地球中心重力・J2・太陽・月による摂動加速度のノルムを高度 (0〜43000 km) に対して片対数プロットする作図スクリプト。慣性系 x 軸方向に宇宙機・月・太陽が並ぶと仮定する。`orbitConst` を呼んで定数を取得する。
- **Note**: 実行可能スクリプト (関数ではない)。`clc`/`clear` を実行し、`const.GM`/`GS`/`GE`/`J2`/`RE` 等の軌道定数に依存する。ヘッダに「途中」とあり作図用の暫定コード。

---

## 姿勢表現・運動学 (`attitude/`)

剛体の姿勢を表現する各種パラメータ（クォータニオン・回転行列 (DCM)・オイラー角・Rodrigues パラメータ・回転ベクトル）とその相互変換、クォータニオン代数（積・共役・逆・平均・補間）、姿勢キネマティクス、TRIAD による姿勢決定をまとめたディレクトリ。

クォータニオンを扱う関数の多くは、第1引数 `scalar` でクォータニオンの成分並びを切り替える共通規約を持つ:
- `scalar == 0`: q = [q0, qv] = [cos(θ/2), e sin(θ/2)] （スカラー先頭）
- `scalar == 4`: q = [qv, q4] = [e sin(θ/2), cos(θ/2)] （スカラー末尾）

積を扱う関数は第2引数 `def` で積の定義を切り替える（`def == 0`: ⊙、`def == 1`: ⊗、後者はベクトル部の外積符号が反転）。ベクトル/クォータニオンは原則 **行方向=サンプル、列方向=成分** の横ベクトル規約（nx3 / nx4）でバッチ処理する。角度単位は rad。`skew`（`math/`）に依存する関数がある。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `dcm1axis` | 単一軸（1/2/3 で指定）まわりの回転行列 |
| `dcm1axisX` | X 軸まわりの回転行列 |
| `dcm1axisY` | Y 軸まわりの回転行列 |
| `dcm1axisZ` | Z 軸まわりの回転行列 |
| `eulerDCM` | 任意の 3 軸オイラー角列 → 回転行列 |
| `zyx2dcm` | ZYX (3-2-1) オイラー角 → 回転行列 |
| `zxz2dcm` | ZXZ (3-1-3) オイラー角 → 回転行列 |
| `zyz2dcm` | ZYZ (3-2-3) オイラー角 → 回転行列 |
| `zyx2q` | ZYX (3-2-1) オイラー角 → クォータニオン |
| `zxz2q` | ZXZ (3-1-3) オイラー角 → クォータニオン |
| `zyz2q` | ZYZ (3-2-3) オイラー角 → クォータニオン |
| `q2zyx` | クォータニオン → ZYX オイラー角 |
| `q2zyz` | クォータニオン → ZYZ オイラー角 |
| `q2dcm` | クォータニオン → 回転行列 |
| `dcm2q` | 回転行列 → クォータニオン（不連続の可能性あり） |
| `dcm2qContinuous` | 回転行列 → 前サンプルと連続なクォータニオン |
| `q1axis` | 単一軸・回転角 → クォータニオン |
| `qAxisAngle` | クォータニオン → オイラー軸・回転角 |
| `q2rotVec` | クォータニオン → 回転ベクトル (θe) |
| `rotVec2q` | 回転ベクトル (θe) → クォータニオン |
| `q2Rodrigues` | クォータニオン → Rodrigues パラメータ |
| `rodrigues2q` | Rodrigues パラメータ → クォータニオン |
| `q2grp` | クォータニオン → 一般化 Rodrigues パラメータ (GRP) |
| `grp2q` | 一般化 Rodrigues パラメータ (GRP) → クォータニオン |
| `qMult` | クォータニオン積 |
| `qMultMat` | クォータニオン積の行列表現 (4x4) |
| `qConj` | クォータニオン共役 |
| `qInv` | クォータニオン逆元 |
| `qErr` | エラークォータニオン（qd を q に一致させる） |
| `qAve` | 重み付きクォータニオン平均 |
| `slerp` | 球面線形補間 (slerp) |
| `qRotation` | クォータニオンによる座標系回転（ベクトルを回転後座標系で表示） |
| `qKine` | クォータニオンキネマティクス（連続時間） |
| `qPropMat` | クォータニオン離散伝播行列 (4x4) |
| `triad` | TRIAD 法による姿勢決定（DCM） |
| `rotPeriod` | 剛体の自由回転運動の周期 |
| `meanAngle` | 重み付き角度平均（円周統計） |
| `generateEulerAngleKinematics` | オイラー角キネマティクス行列 B の記号導出 (スクリプト) |

### 回転行列生成（単一軸・オイラー角）

#### `dcm1axis`
```matlab
R = dcm1axis(axis, phi)
```
単一軸まわりの回転行列。座標変換（フレーム回転）の DCM。
- **入力**: `axis` (1x1, `{1,2,3}`) — 回転軸 (1=x, 2=y, 3=z); `phi` (1x1, rad) — 回転角
- **出力**: `R` (3x3) — 回転行列
- *See also*: `dcm2q`

#### `dcm1axisX` / `dcm1axisY` / `dcm1axisZ`
```matlab
R = dcm1axisX(phi)
R = dcm1axisY(phi)
R = dcm1axisZ(phi)
```
`dcm1axis` の各軸専用版。X/Y/Z 軸まわりの回転行列を返す。
- **入力**: `phi` (1x1, rad) — 回転角
- **出力**: `R` (3x3) — 回転行列
- *See also*: `dcm1axis`

#### `eulerDCM`
```matlab
R = eulerDCM(axis1, axis2, axis3, phi, theta, psi)
```
任意に指定した 3 軸のオイラー角列から回転行列を生成。R = dcm1axis(axis3,psi) · dcm1axis(axis2,theta) · dcm1axis(axis1,phi)。
- **入力**: `axis1,axis2,axis3` (各 1x1, `{1,2,3}`) — 第1/2/3 回転軸; `phi,theta,psi` (各 1x1, rad) — 第1/2/3 回転角
- **出力**: `R` (3x3) — 回転行列
- **Note**: ZYX (3-2-1) が欲しい場合は `eulerDCM(3,2,1, phi,theta,psi)`。

#### `zyx2dcm`
```matlab
R = zyx2dcm(phi, theta, psi)
```
ZYX (3-2-1) オイラー角から回転行列を生成。
- **入力**: `phi` (rad) — Z 軸まわり第1回転; `theta` (rad) — Y 軸まわり第2回転; `psi` (rad) — X 軸まわり第3回転
- **出力**: `R` (3x3) — 回転行列

#### `zxz2dcm`
```matlab
R = zxz2dcm(phi, theta, psi)
```
ZXZ (3-1-3) オイラー角から回転行列を生成。
- **入力**: `phi` (rad) — Z 軸まわり第1回転; `theta` (rad) — X 軸まわり第2回転; `psi` (rad) — Z 軸まわり第3回転
- **出力**: `R` (3x3) — 回転行列
- *See also*: `zyx2dcm`

#### `zyz2dcm`
```matlab
R = zyz2dcm(phi, theta, psi)
```
ZYZ (3-2-3) オイラー角から回転行列を生成。
- **入力**: `phi` (rad) — Z 軸まわり第1回転; `theta` (rad) — Y 軸まわり第2回転; `psi` (rad) — Z 軸まわり第3回転
- **出力**: `R` (3x3) — 回転行列

### オイラー角 ↔ クォータニオン

#### `zyx2q`
```matlab
q = zyx2q(scalar, phi, theta, psi)
```
ZYX (3-2-1) オイラー角からクォータニオンを生成。
- **入力**: `scalar` (1x1, `{0,4}`) — 出力クォータニオンの定義; `phi` (nx1, rad) — Z 軸第1回転; `theta` (nx1, rad) — Y 軸第2回転; `psi` (nx1, rad) — X 軸第3回転
- **出力**: `q` (nx4) — クォータニオン
- *See also*: `qMult`

#### `zxz2q`
```matlab
q = zxz2q(scalar, phi, theta, psi)
```
ZXZ (3-1-3) オイラー角からクォータニオンを生成。
- **入力**: `scalar` (1x1, `{0,4}`); `phi` (nx1, rad) — Z 軸第1回転; `theta` (nx1, rad) — X 軸第2回転; `psi` (nx1, rad) — Z 軸第3回転
- **出力**: `q` (nx4) — クォータニオン
- *See also*: `qMult`

#### `zyz2q`
```matlab
q = zyz2q(scalar, phi, theta, psi)
```
ZYZ (3-2-3) オイラー角からクォータニオンを生成。
- **入力**: `scalar` (1x1, `{0,4}`); `phi` (nx1, rad) — Z 軸第1回転; `theta` (nx1, rad) — Y 軸第2回転; `psi` (nx1, rad) — Z 軸第3回転
- **出力**: `q` (nx4) — クォータニオン
- *See also*: `qMult`, `zyx2q`

#### `q2zyx`
```matlab
euler = q2zyx(scalar, q)
```
クォータニオンから ZYX オイラー角を計算。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `euler` (nx3, rad) — [phi, theta, psi]（phi: z 軸第1回転、theta: y 軸第2回転、psi: x 軸第3回転）

#### `q2zyz`
```matlab
euler = q2zyz(scalar, q)
```
クォータニオンから ZYZ (3-2-3) オイラー角を計算。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `euler` (nx3, rad) — [phi, theta, psi]（phi: z 軸第1回転、theta: y 軸第2回転、psi: z 軸第3回転）

### クォータニオン ↔ 回転行列

#### `q2dcm`
```matlab
DCM = q2dcm(scalar, q)
```
クォータニオンから回転行列を計算（内部で正規化）。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (1x4) — クォータニオン
- **出力**: `DCM` (3x3) — 回転行列
- **Note**: 実装は単一クォータニオン（1 サンプル）を前提としたスカラー要素アクセス。バッチ入力は不可。

#### `dcm2q`
```matlab
q = dcm2q(scalar, R)
```
回転行列からクォータニオンを計算。Markley (2008) の数値安定な手法。
- **入力**: `scalar` (1x1, `{0,4}`); `R` (3x3) — 回転行列
- **出力**: `q` (1x4) — クォータニオン
- **Note**: この方法では時系列でクォータニオンが不連続になりうる。連続性が必要なら `dcm2qContinuous` を使う（コメント中の `dcm2qC.mlx` は現存せず、実体は `dcm2qContinuous.m`）。
- **参照**: Markley, F. L. (2008). "Unit Quaternion from Rotation Matrix," JGCD 31, 440-442.

#### `dcm2qContinuous`
```matlab
q = dcm2qContinuous(scalar, qK, dcm)
```
前サンプルのクォータニオンに符号連続なクォータニオンを回転行列から計算。Wu (2019) の最適連続化。
- **入力**: `scalar` (1x1, `{0,4}`); `qK` (1x4) — 前サンプル q_{k-1}; `dcm` (3x3) — 回転行列
- **出力**: `q` (1x4) — 連続化クォータニオン
- **参照**: Wu, J. (2019). "Optimal Continuous Unit Quaternions from Rotation Matrices," JGCD 42(4), 919-922.
- *See also*: `dcm2q`

### クォータニオン ↔ 軸角・回転ベクトル・Rodrigues

#### `q1axis`
```matlab
q = q1axis(scalar, axis, theta)
```
単一軸と回転角からクォータニオンを生成（軸は内部で正規化）。
- **入力**: `scalar` (1x1, `{0,4}`); `axis` (nx3) — 回転軸（単位ベクトル、非正規でも可）; `theta` (nx1, rad) — 回転角
- **出力**: `q` (nx4) — クォータニオン
- *See also*: `dcm1axis`

#### `qAxisAngle`
```matlab
[eAxis, eAngle] = qAxisAngle(scalar, q)
```
クォータニオンからオイラー軸と回転角を抽出。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `eAxis` (nx3) — オイラー軸 e; `eAngle` (nx1, rad) — 回転角 θ ∈ [0, 2π)（mod 2π で丸め）
- *See also*: `qConj`, `qInv`

#### `q2rotVec`
```matlab
rotVec = q2rotVec(scalar, q)
```
クォータニオンから回転ベクトル θe を計算（内部で `qAxisAngle` を利用）。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `rotVec` (nx3) — 回転ベクトル θe

#### `rotVec2q`
```matlab
q = rotVec2q(scalar, rv)
```
回転ベクトル θe からクォータニオンを計算。
- **入力**: `scalar` (1x1, `{0,4}`); `rv` (nx3) — 回転ベクトル θe
- **出力**: `q` (nx4) — クォータニオン

#### `q2Rodrigues`
```matlab
Rod = q2Rodrigues(scalar, q)
```
クォータニオンから（古典的）Rodrigues パラメータを計算。内部で `q2grp` を f=1, a=0 で呼ぶ。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `Rod` (nx3) — Rodrigues パラメータ

#### `rodrigues2q`
```matlab
q = rodrigues2q(scalar, rod)
```
Rodrigues パラメータからクォータニオンを計算。内部で `grp2q` を f=1, a=0 で呼ぶ。
- **入力**: `scalar` (1x1, `{0,4}`); `rod` (nx3) — Rodrigues パラメータ
- **出力**: `q` (nx4) — クォータニオン

#### `q2grp`
```matlab
p = q2grp(scalar, f, a, q)
```
クォータニオンから一般化 Rodrigues パラメータ (GRP) を計算。
- **入力**: `scalar` (1x1, `{0,4}`); `f` (1x1) — スケーリングパラメータ; `a` (1x1) — GRP パラメータ; `q` (nx4) — クォータニオン
- **出力**: `p` (nx3) — 一般化 Rodrigues パラメータ

#### `grp2q`
```matlab
q = grp2q(scalar, f, a, p)
```
一般化 Rodrigues パラメータ (GRP) からクォータニオンを計算。
- **入力**: `scalar` (1x1, `{0,4}`); `f` (1x1) — スケーリングパラメータ; `a` (1x1) — GRP パラメータ; `p` (nx3) — 一般化 Rodrigues パラメータ
- **出力**: `q` (nx4) — クォータニオン

### クォータニオン代数

#### `qMult`
```matlab
output = qMult(scalar, def, q, p)
```
クォータニオン積。バッチ（要素ごと）演算に対応。
- **入力**: `scalar` (1x1, `{0,4}`); `def` (1x1, `{0,1}`) — 積の定義（0: ⊙、1: ⊗、ベクトル部外積符号が反転）; `q` (nx4); `p` (nx4)
- **出力**: `output` (nx4) — 積クォータニオン
- *See also*: `qConj`, `qInv`

#### `qMultMat`
```matlab
qMat = qMultMat(scalar, def, q)
```
クォータニオン積を行列 · ベクトル形式にするための 4x4 行列表現。
- **入力**: `scalar` (1x1, `{0,4}`); `def` (1x1, `{0,1}`) — 積の定義; `q` (1x4) — クォータニオン（単一）
- **出力**: `qMat` (4x4) — 積行列
- **Note**: 単一クォータニオン前提（スカラー要素アクセス）。バッチ入力は不可。
- *See also*: `qConj`, `qInv`

#### `qConj`
```matlab
qC = qConj(scalar, q)
```
クォータニオン共役（ベクトル部の符号反転）。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `qC` (nx4) — 共役
- *See also*: `qInv`

#### `qInv`
```matlab
qInv = qInv(scalar, q)
```
クォータニオン逆元 q^{-1} = q* / ‖q‖^2（正規化を含む）。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — クォータニオン
- **出力**: `qInv` (nx4) — 逆元
- *See also*: `qConj`

#### `qErr`
```matlab
qe = qErr(scalar, q, qd)
```
qd を回転させて q に一致させるためのエラークォータニオン。qe = q ⊗ qd^{-1}（内部で `qInv` + `qMult(def=1)`）。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (nx4) — 目標（真値）; `qd` (nx4) — desired / 推定クォータニオン
- **出力**: `qe` (nx4) — エラークォータニオン

#### `qAve`
```matlab
qAveraged = qAve(q, w)
```
重み付きクォータニオン平均。M = Σ w_i q_i q_i^T の最大固有ベクトル（Markley 2007）。出力の定義は入力と整合。
- **入力**: `q` (nx4) — クォータニオン群; `w` (nx1) — スカラー重み（省略時は一様重み 1/n）
- **出力**: `qAveraged` (1x4) — 平均クォータニオン
- **Note**: `scalar` 引数を取らない（入力定義がそのまま維持される）。
- **参照**: Markley et al. (2007). "Averaging Quaternions," JGCD 30(4), 1193-1197.

#### `slerp`
```matlab
qt = slerp(t, scalar, q1, q2)
```
2 つのクォータニオン間の球面線形補間。
- **入力**: `t` (正規化時刻, 0 ≤ t ≤ 1、ベクトル可); `scalar` (1x1, `{0,4}`); `q1` (1x4) — 始点; `q2` (1x4) — 終点
- **出力**: `qt` — 補間クォータニオン（`t` の要素数 x 4）
- *See also*: `sclerp`（`dualQuaternions/`、デュアルクォータニオン版）

#### `qRotation`
```matlab
rb = qRotation(scalar, r, q)
```
ベクトルを回転後の座標系で表示（ベクトル自体を回転させるのではなくフレーム回転）。r_b = q^{-1} ⊙ r ⊙ q。
- **入力**: `scalar` (1x1, `{0,4}`); `r` (nx3) — ベクトル; `q` (nx4) — クォータニオン
- **出力**: `rb` (nx3) — 回転後座標系で表したベクトル
- *See also*: `qInv`

### キネマティクス

#### `qKine`
```matlab
qKine = qKine(scalar, q, w)
```
連続時間クォータニオンキネマティクス q̇ = (1/2) ω ⊗ q を計算。
- **入力**: `scalar` (1x1, `{0,4}`); `q` (1x4) — クォータニオン（内部で列ベクトル化）; `w` (1x3, rad/s) — 角速度
- **出力**: `qKine` (1x4) — クォータニオン時間微分 q̇
- **Note**: 単一状態前提。`scalar` により運動学行列が切り替わる。

#### `qPropMat`
```matlab
qProp = qPropMat(scalar, dt, w)
```
一定角速度を仮定した離散クォータニオン伝播行列（q_{k+1} = qProp · q_k に用いる）。
- **入力**: `scalar` (1x1, `{0,4}`); `dt` (1x1, s) — 時間刻み; `w` (1x3, rad/s) — 角速度
- **出力**: `qProp` (4x4) — 離散伝播行列
- **Note**: 出力はクォータニオンではなく 4x4 行列。`skew`（`math/`）に依存。
- *See also*: `qInv`

#### `generateEulerAngleKinematics`
```matlab
% スクリプト（引数・戻り値なし）
```
オイラー角キネマティクス φ̇θ̇ψ̇ = B(θ,ψ) ω の変換行列 B を Symbolic Math Toolbox で記号導出するスクリプト。冒頭 `first/second/third` に回転軸 (1/2/3) を設定して実行し、簡約化した B とその LaTeX を出力する。
- **Note**: 関数ではなくスクリプト（`clc; clear` を含む）。実際の出力は記号行列 `B`（および `latex(B)`）で、`dxdt` という変数は生成しない。Symbolic Math Toolbox が必要。

### 姿勢決定・その他

#### `triad`
```matlab
R = triad(v1_i, v2_i, w1_b, w2_b)
```
2 組の参照ベクトルから TRIAD 法で慣性系→機体固定系の回転行列を計算（各ベクトルは内部で正規化）。
- **入力**: `v1_i, v2_i` (各 nx3) — 慣性系で表した参照ベクトル; `w1_b, w2_b` (各 nx3) — 機体固定系で表した参照ベクトル
- **出力**: `R` (3x3xn) — 慣性系→機体固定系の回転行列（サンプルごとに 3 枚目次元へ格納）
- *See also*: `q2dcm`, `dcm2q`

#### `rotPeriod`
```matlab
T = rotPeriod(MOI, w)
```
剛体の自由回転運動（トルクフリー）の周期を第一種完全楕円積分から計算。
- **入力**: `MOI` (3x3, kg·m^2) — 主軸慣性モーメント（対角成分 Jx,Jy,Jz を使用）; `w` (1x3, rad/s) — 角速度
- **出力**: `T` (s) — 自由回転周期

#### `meanAngle`
```matlab
thetaAve = meanAngle(thetaArray, w)
```
角度の重み付き平均（円周統計）。cos/sin を重み付き平均して atan2 で復元し、周期性 (±π ラップ) を正しく扱う。
- **入力**: `thetaArray` (ベクトル, rad) — 角度群; `w` — 重み（内部で正規化）
- **出力**: `thetaAve` (rad) — 平均角度 ∈ (-π, π]
- **Note**: ドキュメントコメントを持たない裸関数。単純な算術平均ではなく円周平均である点に注意（境界をまたぐ角度で結果が変わる）。

---

## デュアルクォータニオン (`dualQuaternions/`)

剛体の姿勢と位置をまとめて表現するデュアルクォータニオン（dual quaternion）の生成・変換・演算・補間を提供する。すべての関数が横ベクトル規約（行方向=時系列/サンプル、列方向=状態成分）に従い、real part と dual part を連結した `nx8` 行列で dual quaternion を扱う（`dq = [real(1:4), dual(5:8)]`）。多くの関数は先頭引数 `scalar` でクォータニオンのスカラ部位置の定義を切り替える（`scalar==0` は先頭 q0、`scalar==4` は末尾 q4）。回転・並進の演算は `attitude/` の `qMult` / `qConj` / `qInv` / `qMultMat` に依存する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `pos2dq` | 位置ベクトルとクォータニオンを dual quaternion へ変換 |
| `dq2pos` | dual quaternion を位置ベクトルとクォータニオンへ逆変換 |
| `dqMult` | 2 つの dual quaternion の積 |
| `dqMultMat` | dual quaternion 積を表す 8x8 行列を構成 |
| `dqConj` | 共役 dual quaternion |
| `dqInv` | 逆 dual quaternion |
| `cay` | ベクトルからクォータニオンへのケイリー変換 |
| `cayInv` | クォータニオンからベクトルへの逆ケイリー変換 |
| `dqCay` | dual ベクトルの dual ケイリー変換（dual quaternion を生成） |
| `sclerp` | dual quaternion 間のスクリュー線形補間 |
| `ctrlDq` | C1 連続 sclerp のための制御 dual quaternion を生成 |

### Core（変換・演算）

#### `pos2dq`
```matlab
dq = pos2dq(inertial, scalar, r, q)
```
位置ベクトル `r` と姿勢クォータニオン `q` を dual quaternion に変換する。
- **入力**: `inertial` (1x1, {0,1}) — 位置ベクトルの表現フレーム。0: body-fixed frame、1: inertial frame; `scalar` (1x1, {0,4}) — スカラ部の定義; `r` (nx3) — 位置ベクトル; `q` (nx4) — クォータニオン
- **出力**: `dq` (nx8) — dual quaternion `[real, dual]`
- *See also*: `qMult`

#### `dq2pos`
```matlab
[r, q] = dq2pos(inertial, scalar, dq)
```
dual quaternion を位置ベクトルとクォータニオンに逆変換する（`pos2dq` の逆）。
- **入力**: `inertial` (1x1, {0,1}) — 出力位置ベクトルの表現フレーム。1: inertial frame、0: body-fixed frame; `scalar` (1x1, {0,4}) — スカラ部の定義; `dq` (nx8) — dual quaternion `[real, dual]`
- **出力**: `r` (nx3) — 位置ベクトル; `q` (nx4) — クォータニオン（real part そのもの）
- *See also*: `dqMult`, `pos2dq`

#### `dqMult`
```matlab
dq = dqMult(scalar, def, dq, dp)
```
2 つの dual quaternion の積を計算する。real part 同士の積と、real×dual の交差項からなる dual part を返す。
- **入力**: `scalar` (1x1, {0,4}) — スカラ部の定義; `def` (1x1, {0,1}) — クォータニオン積の定義。0: 右手系（×+）、1: 反対の符号規約（×−）; `dq` (nx8) — dual quaternion `[real, dual]`; `dp` (nx8) — dual quaternion `[real, dual]`
- **出力**: `dq` (nx8) — 積の dual quaternion（入力 `dq` を上書きする名前）
- *See also*: `qMult`

#### `dqMultMat`
```matlab
dqM = dqMultMat(scalar, def, dq)
```
左からの dual quaternion 積 `dq ⊙ (·)` を線形写像として表す 8x8 行列を構成する。ブロック下三角構造で、対角に real part、左下に dual part の積行列を配置する。
- **入力**: `scalar` (1x1, {0,4}) — スカラ部の定義; `def` (1x1, {0,1}) — クォータニオン積の定義; `dq` (1x8) — dual quaternion `[real, dual]`
- **出力**: `dqM` (8x8) — dual quaternion 積行列
- **Note**: 内部で `qMultMat` を呼ぶため入力は 1 つの dual quaternion（1x8）を想定する（バッチ入力は非対応）。
- *See also*: `qMult`, `qMultMat`, `dqMult`

#### `dqConj`
```matlab
invDq = dqConj(scalar, dq)
```
共役 dual quaternion を返す。real/dual それぞれのベクトル部の符号を反転する。
- **入力**: `scalar` (1x1, {0,4}) — スカラ部の定義; `dq` (nx8) — dual quaternion `[real, dual]`
- **出力**: `invDq` (nx8) — 共役 dual quaternion `[real, dual]`
- **Note**: 出力変数名が `invDq` だが返すのは共役（conjugate）であって逆元ではない。逆元は `dqInv` を使う。
- *See also*: `qMult`

#### `dqInv`
```matlab
invDq = dqInv(scalar, dq)
```
逆 dual quaternion を返す。real part は `qInv`、dual part は `-qr^{-1} ⊙ qd ⊙ qr^{-1}` で構成する。
- **入力**: `scalar` (1x1, {0,4}) — スカラ部の定義; `dq` (nx8) — dual quaternion `[real, dual]`
- **出力**: `invDq` (nx8) — 逆 dual quaternion `[real, dual]`
- *See also*: `qMult`, `dqConj`

### Cayley 変換

#### `cay`
```matlab
q = cay(scalar, u)
```
ケイリー変換 cay(u) = (1+u) ⊙ (1−u)^{-1}。任意ベクトルからクォータニオンを生成する。
- **入力**: `scalar` (1x1, {0,4}) — 出力クォータニオンのスカラ部定義; `u` (nx3) — 任意のベクトル
- **出力**: `q` (nx4) — クォータニオン
- **Note**: 内部計算はスカラ部末尾（q4）の規約で行い、`scalar==0` 指定時は最後に先頭スカラ部形式へ並べ替える。
- *See also*: `qMult`

#### `cayInv`
```matlab
u = cayInv(scalar, q)
```
逆ケイリー変換 cay^{-1}(q) = (q−1) ⊙ (q+1)^{-1}。クォータニオンから対応するベクトルを取り出す。
- **入力**: `scalar` (1x1, {0,4}) — 入力クォータニオンのスカラ部定義; `q` (nx4) — クォータニオン
- **出力**: `u` (nx3) — ベクトル
- **Note**: 定義式の「1」は単位クォータニオン（スカラ部=1、ベクトル部=0）。内部で末尾スカラ部形式に統一して計算する。
- *See also*: `qMult`

#### `dqCay`
```matlab
dq = dqCay(scalar, dv)
```
dual ベクトルの dual ケイリー変換 cay(ũ) = (1+ũ) ⊙ (1−ũ)^{-1}。real part は `cay`、dual part は解析式で計算し、dual quaternion を生成する。
- **入力**: `scalar` (1x1, {0,4}) — スカラ部の定義; `dv` (nx8) — 任意の dual ベクトル（スカラ部は 0 を想定）`[real, dual]`
- **出力**: `dq` (nx8) — dual quaternion
- *See also*: `cay`

### 補間

#### `sclerp`
```matlab
dqt = sclerp(t, scalar, dq1, dq2)
```
2 つの dual quaternion 間のスクリュー線形補間（ScLERP）。ねじ運動軸まわりに回転と並進を同時に補間する。
- **入力**: `t` (nx1) — 正規化時刻（0 < t ≤ 1、列ベクトル化される）; `scalar` (1x1, {0,4}) — スカラ部の定義; `dq1` (1x8) — 始点 dual quaternion; `dq2` (1x8) — 終点 dual quaternion
- **出力**: `dqt` (nx8) — 補間された dual quaternion（`t` のサンプル数 n 行）
- **Note**: 回転がほぼゼロ（sin(θ/2) が eps 未満）の場合は始点 `dq1` を全時刻に複製して返す。`dq1`/`dq2` は 1 サンプル（1x8）を前提とする。
- *See also*: （明示なし）

#### `ctrlDq`
```matlab
[dqa, dqb] = ctrlDq(dtp, dq1, dq2, w1, v1, w2, v2)
```
C1 連続な sclerp 曲線（Allmendinger et al. 2018）のための 2 つの制御 dual quaternion を、端点の角速度・並進速度から構成する。
- **入力**: `dtp` (1x1, s) — 時間ステップ; `dq1` (1x8) — 補間始点 dual quaternion; `dq2` (1x8) — 補間終点 dual quaternion; `w1` (1x3) — 始点の角速度; `v1` (1x3) — 始点の並進速度; `w2` (1x3) — 終点の角速度; `v2` (1x3) — 終点の並進速度
- **出力**: `dqa` (1x8) — 始点側制御 dual quaternion; `dqb` (1x8) — 終点側制御 dual quaternion
- **Note**: 実際のシグネチャは `dq1, dq2` の後に `w1, v1, w2, v2` の順（端点ごとに角速度・並進速度が交互）で並ぶ。ドキュメントコメントの記述順（w1,w2 / v1,v2）とは異なるので呼び出し順に注意。単位は `dtp` と速度で整合させること。
- *See also*: `sclerp`, `qMult`, `dqConj`

---

## 幾何学的積分 (`geometricIntegration/`)

剛体姿勢のクォータニオンを Lie 群構造を保ったまま時間発展させる幾何学的数値積分（Crouch–Grossman 法）と、その検証に使う軸対称剛体の解析解・剛体運動方程式・Runge–Kutta / Crouch–Grossman の Butcher 係数表をまとめたディレクトリです。時系列出力は **行方向=時刻ステップ、列方向=状態成分** の横ベクトル規約（`qOut` は N×4、`wOut` は N×3）で返します。クォータニオンの成分順序は `scalar` 引数で切り替え、`scalar==4` はスカラー部を末尾に置く `[e·sin(θ/2), cos(θ/2)]` 規約、`scalar==0` はスカラー部を先頭に置く `[cos(θ/2), e·sin(θ/2)]` 規約です。角速度・慣性テンソル（`MOI`、3×3）に依存し、内部で `attitude/` の `qMult` / `qMultMat` を利用します。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `qGI` | Crouch–Grossman 法によるクォータニオンの幾何学的積分（トルクフリー剛体） |
| `axiQsol` | 軸対称剛体のトルクフリー姿勢運動の解析解（Andrle & Crassidis 2013） |
| `qExp` | クォータニオン指数写像の閉形式（4×4 行列） |
| `eulerEom` | 剛体の姿勢運動方程式（ode ソルバ用の右辺、トルクフリー） |
| `butcherTable` | RK 法 / Crouch–Grossman 法の Butcher 係数表 (a, b, c) |

### Core

#### `qGI`
```matlab
[qOut, wOut] = qGI(scalar, tspan, qIni, wIni, nGI, MOI)
```
クォータニオン姿勢のトルクフリー剛体運動を Crouch–Grossman 法（Lie 群上の陽的積分）で時間発展させる mixed-scheme 実装。各ステージで角速度を Runge–Kutta 的に更新しつつ、クォータニオンは `expm` による行列指数で更新する。
- **入力**: `scalar` (1,1) — クォータニオン成分順序の指定（`scalar==1` のとき入力 `qIni` はスカラー先頭とみなし内部でスカラー末尾へ並べ替える。それ以外はそのまま使う）; `tspan` (1,N) — 時刻列（等間隔である必要はなく `diff` でステップ幅を取る）; `qIni` (1,4) — 初期クォータニオン; `wIni` (1,3, rad/s) — 初期角速度; `nGI` (1,1) — 積分次数、3 または 4（それぞれ 3 ステージ / 5 ステージの CG 法に対応）; `MOI` (3,3, kg·m^2) — 慣性テンソル
- **出力**: `qOut` (N,4) — 各時刻のクォータニオン（スカラー末尾規約）; `wOut` (N,3, rad/s) — 各時刻の角速度
- **Note**: 係数は `butcherTable(nGI,'CG')` から取得するため `nGI` は 3 か 4 のみ。トルクは 0（外力なし）で固定されている。クォータニオン更新は既定で `expm`、コメントアウトされた `qExp` による閉形式が代替として用意されている。
- *See also*: `butcherTable`, `qExp`, `qMultMat`

### Verification / Reference solutions

#### `axiQsol`
```matlab
[qOut, wOut] = axiQsol(t_, qIni, wIni, MOI)
```
軸対称剛体のトルクフリー姿勢運動の解析解（Andrle & Crassidis, 2013）。角運動量ベクトルまわりの規則歳差から角速度とクォータニオンの厳密解を与え、`qGI` 等の数値積分器の検証基準として使う。
- **入力**: `t_` (N,1 に整形, s) — 評価時刻ベクトル（内部で列ベクトル化）; `qIni` (1,4) — 初期クォータニオン（スカラー末尾規約で `qMult(4,1,...)` を通して合成）; `wIni` (1,3 または 3,1, rad/s) — 初期角速度; `MOI` (3,3, kg·m^2) — 軸対称慣性テンソル（`MOI(1,1)=MOI(2,2)` を横軸、`MOI(3,3)` を対称軸として使用）
- **出力**: `qOut` (N,4) — 各時刻のクォータニオン; `wOut` (N,3, rad/s) — 各時刻の角速度
- **Note**: 慣性テンソルが軸対称（`MOI(1,1)==MOI(2,2)`）であることが前提。歳差角速度 `wn` は `wIni(3)*(MOI(1,1)-MOI(3,3))/MOI(1,1)` で計算される。**実装は `qMult_` という存在しない関数を呼んでおり、そのままでは実行時にエラーになる（正しくは `qMult`）。**

#### `qExp`
```matlab
qMat = qExp(scalar, dt, w)
```
角速度が一定という仮定のもとでのクォータニオン運動学の指数写像を閉形式で与える 4×4 行列。`eye(4)*cos(dt*|w|/2) + qMultMat*sin(dt*|w|/2)/|w|` の形で、`qGI` のクォータニオン更新の代替（`expm` 相当の閉形式）として使える。
- **入力**: `scalar` (1,1) — クォータニオン成分順序の指定（`qMultMat` に渡される）; `dt` (1,1, s) — 時間刻み; `w` (3,1 に整形, rad/s) — 角速度（内部で列ベクトル化し、ノルムで正規化）
- **出力**: `qMat` (4,4) — クォータニオン更新行列。`qNext = qMat * qNow` の形で作用する
- **Note**: ドキュメントコメントがテンプレート雛形のまま（「この関数の簡単な概要です／詳細な説明です」）で、実際の指数写像の説明になっていない。`w` がゼロベクトルの場合 `wNorm=0` で 0 除算が生じる。

### Rigid-body dynamics

#### `eulerEom`
```matlab
dxdt = eulerEom(t, x, MOI)
```
剛体姿勢の運動方程式の右辺（`ode45` 等のソルバに渡す関数ハンドル用）。クォータニオン運動学と Euler の回転方程式を組み合わせた 7 状態の時間微分を返す。
- **入力**: `t` (1,1, s) — 時刻（自律系のため未使用だがソルバ規約で受ける）; `x` (7,1) — 状態ベクトル `[q(4); w(3)]`（`q` はスカラー末尾規約、`w` は rad/s）; `MOI` (3,3, kg·m^2) — 慣性テンソル
- **出力**: `dxdt` (7,1) — 状態微分 `[dq/dt(4); dw/dt(3)]`
- **Note**: 制御トルク `u` はコード内で `zeros(3,1)` に固定（トルクフリー）で、他の候補式はコメントアウトされている。ヘッダコメントの Input/Output 欄が空で、`used in eulerMain.m` とあるが `eulerMain.m` はライブラリに存在しない。

### Coefficient tables

#### `butcherTable`
```matlab
[a, b, c] = butcherTable(n, method)
```
幾何学的積分に用いる Butcher 係数表を返す。Runge–Kutta 法（`'RK'`）と Crouch–Grossman 法（`'CG'`）に対応し、次数 `n` は 3 と 4 をサポートする。
- **入力**: `n` (1,1) — 積分法の次数（ステージ数ではない）。3 または 4 のみ有効、それ以外は `error`; `method` (char) — `'RK'`（Runge–Kutta）または `'CG'`（Crouch–Grossman）。それ以外は `error`
- **出力**: `a` (s,s) — 段係数行列（RK3/CG3 は 3×3、RK4 は 4×4、CG4 は 5×5）; `b` (s,1) — 重み係数ベクトル; `c` (s,1) — ノード係数ベクトル
- **Note**: `n` は次数でありステージ数と一致しない（CG4 は 4 次だが 5 ステージ、係数配列は 5×5）。ヘッダの出力見出しが `## otuputs` と綴り誤り。`See also dqrk, dqgi` の参照先 `dqrk` / `dqgi` はライブラリに存在しない。参考文献欄は "to be added" のまま未記載。

---

## 宇宙機形状モデル (OBJ) (`object/`)

Wavefront `.obj` / `.mtl` 形式で定義した宇宙機のポリゴンメッシュを読み込み、面ごとの幾何量（法線・面積・重心・ローカルフレーム）や表面光学係数を計算し、自己遮蔽 (self-shadowing) 判定や 3D 可視化を行う関数群。中心となるデータは構造体 `sat` で、面 (facet) を行方向にとる横ベクトル規約（各フィールドは `N×3` 行列または `N×1` ベクトル、`N` = 面数）で統一されている。読み込みは `readSC` を入口とし、内部で `readObj` → `calcAreaObj` → `calcLocalFrame` を呼んでパイプライン的に `sat` を構築する。

**`sat` 構造体の主なフィールド**（パイプラインが実際に生成する名前）:
- `sat.vertices` (N×3, m) — 頂点座標
- `sat.faces` (N×3 または N×4) — 面を構成する頂点インデックス（三角・四角混在時は四角の欠損列を `NaN` で埋める）
- `sat.normal` (N×3) — 面の外向き単位法線ベクトル
- `sat.area` (N×1, m^2) — 面積
- `sat.pos` (N×3, m) — 面の重心
- `sat.Ca` / `sat.Cd` / `sat.Cs` (N×1) — 吸収率・拡散反射率・鏡面反射率（MTL の Kd の R/G/B を割り当て）
- `sat.uu` / `sat.uv` / `sat.qlb` — ローカルフレーム軸とボディ→ローカルのクォータニオン
- `sat.sunlitFlag` (N×1) — 自己遮蔽フラグ（1: 照射、0: 影）

**注意（ドキュメントと実装の不一致）**: いくつかの関数のヘッダコメントは法線フィールドを `sat.n` と記載しているが、実装が実際に読み書きするのは `sat.normal` である（`readObj`/`readSC`/`showSC`/`selfShadow`/`calcLocalFrame` はすべて `sat.normal`）。本文では実装に合わせて `sat.normal` を正とする。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `readSC` | `.obj`/`.mtl` を読み込み、幾何量・光学係数・物理量まで揃えた `sat` を構築する（入口関数） |
| `readObj` | Wavefront `.obj` をパースし頂点・面・法線・マテリアルインデックスを取得 |
| `calcAreaObj` | 各面の面積と重心を計算（三角・四角メッシュ混在対応） |
| `calcNormalObj` | 各面の外向き単位法線ベクトルを計算 |
| `calcLocalFrame` | 各面のローカルフレーム軸とボディ→ローカルのクォータニオンを計算 |
| `selfShadow` | 太陽方向に対し全面の自己遮蔽を判定し `sat.sunlitFlag` を設定 |
| `calcRayIntersect` | 1 本の光線と 1 枚の三角面の交差判定（self-shadowing の下請け） |
| `showSC` | 読み込んだ宇宙機形状を patch で 3D 表示（法線・ローカル x 軸の矢印も任意描画） |

### Core（読み込み・構築）

#### `readSC`
```matlab
sat = readSC(satName)
```
衛星形状ファイルと表面特性を読み込み、幾何量・光学係数・既定の物理量を格納した `sat` 構造体を返す入口関数。内部で `readObj`・`calcAreaObj`・`calcLocalFrame` を呼び、同じディレクトリに `.mtl` があれば各マテリアルの Kd (R,G,B) を面ごとの (Ca, Cd, Cs) に割り当てる。
- **入力**: `satName` (char) — `.obj` ファイルのパス
- **出力**: `sat` (struct) — `vertices`/`faces`/`normal`/`area`/`pos`/`Ca`/`Cd`/`Cs`/`uu`/`uv`/`qlb` に加え、既定値として `F0`, `kappa`, `nu`, `nv`, `mCT`（BRDF モデル用）、`MOI`, `m`（慣性・質量）、`sunlitFlag`, `force`, `torque` などを付与
- **Note**: MTL の RGB を吸収・拡散・鏡面係数として流用する規約（R→Ca, G→Cd, B→Cs）。MTL が無い/見つからない場合は既定値（Cd=Cs=0.5 等）のまま。読み込み時に頂点数・面数を `disp` 表示する。
- *See also*: `showSC`, `readObj`, `calcAreaObj`, `calcLocalFrame`

#### `readObj`
```matlab
obj = readObj(fname)
```
Wavefront `.obj` ファイルをパースし、頂点・面インデックス・面法線・マテリアル情報を取り出す。三角面・四角面の混在に対応し、四角に満たない面の欠損頂点列は `NaN` で埋める。
- **入力**: `fname` (char) — `.obj` ファイルのフルパス
- **出力**: `obj` (struct) — `vertices` (N×3), `faces` (N×3 または N×4), `normal` (N×3), `mtlFileName`, `materialNames`, `materialIndices`
- **Note**: B. Abayowa (Tec^Edge, 2007) の `readObj` を Yoshimura が改変したもの。`usemtl` を追跡して面ごとの材質インデックスを付与する。テクスチャ座標 `vt` はパースするが出力には含めない（`obj.vt` はコメントアウト）。
- *See also*: `showSC`, `readSC`

#### `calcAreaObj`
```matlab
[area, pos] = calcAreaObj(sat)
```
各面の面積と重心（頂点平均座標）を計算する。三角メッシュは外積ノルムの半分、四角メッシュは外積ノルムをそのまま面積とする。三角・四角の混在を各行ごとに判定する。
- **入力**: `sat` (struct) — `sat.vertices` (N×3, m), `sat.faces` (N×3 または N×4)
- **出力**: `area` (N×1, m^2) — 面積; `pos` (N×3, m) — 面の重心
- **Note**: 三角/四角判定は `sat.faces(i,4)` が `NaN` かどうかで行う（4 列目が無い、または `NaN` なら三角）。四角面の面積は対角で分割せず単一外積のノルムで近似している点に注意（平行四辺形近似）。
- *See also*: `showSC`, `readSC`

#### `calcNormalObj`
```matlab
n = calcNormalObj(sat)
```
各面の外向き単位法線ベクトルを、面の頂点 1→2 と 1→3 のベクトルの外積から計算する。長さが極小（≤1e-8）の面はゼロ割回避の微小項を加える。
- **入力**: `sat` (struct) — `sat.vertices` (N×3, m), `sat.faces` (N×3)
- **出力**: `n` (N×3) — 面の外向き単位法線ベクトル
- **Note**: ヘッダコメントは `sat.n` を入力フィールドとして列挙しているが、実際に読むのは `vertices` と `faces` のみで、法線は本関数の出力。呼び出し側は結果を `sat.normal` に格納する運用（`readObj` は独自に法線を読むため、本関数は再計算用途）。
- *See also*: `showSC`, `readSC`, `calcAreaObj`

#### `calcLocalFrame`
```matlab
[uu, uv, qlb] = calcLocalFrame(sat)
```
各面のローカルフレーム（x 軸 = 頂点 1→2 方向、z 軸 = 面法線）を求め、TRIAD 法で方向余弦行列を作り、ボディ固定座標系からローカル座標系へのクォータニオンを計算する。
- **入力**: `sat` (struct) — `sat.faces`, `sat.vertices` (N×3, m), `sat.normal` (N×3)
- **出力**: `uu` (N×3) — ローカル x 軸（ボディ系表現）; `uv` (N×3) — ローカル y 軸 (= normal × uu); `qlb` (N×4) — ボディ→ローカルのクォータニオン
- **Note**: 実装は `sat.normal` を参照するが、ヘッダは法線を `sat.n` と誤記している。TRIAD と `dcm2q`（スカラ部先頭指定 `dcm2q(4, ...)`）は `attitude/` の関数に依存する。
- *See also*: `showSC`, `readSC`, `calcAreaObj`

### Self-shadowing（自己遮蔽）

#### `selfShadow`
```matlab
sat = selfShadow(sat, sun) %#codegen
```
太陽方向ベクトルに対し、太陽光を受ける全面について他面による自己遮蔽を判定し、`sat.sunlitFlag`（1: 照射、0: 影）を設定する。四角面は 2 つの三角に分けて 2 回判定する。`parfor` で並列化（Parallel Computing Toolbox 非導入時は `for` に置換）。
- **入力**: `sat` (struct) — `sat.vertices`, `sat.faces`, `sat.normal`, `sat.pos`; `sun` (1×3 または 3×1) — ボディ固定系での太陽方向ベクトル（内部で正規化）
- **出力**: `sat` (struct) — `sat.sunlitFlag` (N×1) を追加・更新
- **Note**: 太陽側を向く面 (`sat.normal * sun > 0`) のみ遮蔽計算し、他は既定で影扱い。四角面判定は先頭面の 4 列目 `faces(1,4)` の `NaN` 有無でメッシュ全体の三角/四角を分岐しており、面ごとに三角・四角が混在するモデルでは意図通り動かない可能性がある（三角・四角混在時の要注意点）。`See also selfShadow` はヘッダの自己参照になっており、意図は `calcRayIntersect` と思われる。
- *See also*: `calcRayIntersect`（ヘッダ記載は自己参照）

#### `calcRayIntersect`
```matlab
flag = calcRayIntersect(sun, nJ, vertJ, vertI) %#codegen
```
太陽方向の光線が、面 i の重心 `vertI` から見て面 j（頂点 `vertJ`、法線 `nJ`）に遮られるかを判定する下請け関数。交点を求めて三角内外を外積の符号一致で判定する。
- **入力**: `sun` (1×3, 正規化済み想定) — 太陽方向; `nJ` (1×3) — 面 j の法線; `vertJ` (3×3) — 面 j の 3 頂点座標; `vertI` (3×1 または 1×3, m) — 面 i の評価点（重心）
- **出力**: `flag` — 1: 遮られない（照射）、0: 面 j が影を作る
- **Note**: ヘッダコメントは `satName` と `sat` 構造体一式（vertices/faces/area/pos/Ca/Cd/Cs/n/shadowFlag）を入力として列挙しているが、実際のシグネチャは `sun, nJ, vertJ, vertI` の 4 引数でありヘッダと全く対応していない（`arguments` ブロックが正しい仕様）。`selfShadow` からのみ呼ばれる。
- *See also*: `showSC`

### 可視化

#### `showSC`
```matlab
showSC(sat)
showSC(sat, Normal="on")
```
読み込んだ宇宙機形状を `patch` で 3D 表示する。オプションで各面の法線ベクトルとローカル x 軸を `quiver3` の矢印で重ね描きする。
- **入力**: `sat` (struct) — `sat.faces`, `sat.vertices`；`Normal="on"` 指定時は `sat.pos`, `sat.normal`, `sat.uu` も使用; `options.Normal` (char, `'on'`/`'off'`, 既定 `'off'`) — 法線・ローカル軸の矢印描画
- **出力**: なし（figure を生成）
- **Note**: 法線矢印描画には `sat.normal` と `sat.uu` を用いる。ヘッダは `sat.n` および `sat.Ca/Cd/Cs` を列挙するが、後者は本関数では未使用。
- *See also*: `readSC`

---

## 太陽輻射圧 (SRP) (`srp/`)

この節は、衛星に働く太陽輻射圧 (SRP) の力・トルクを反射モデル別に計算する関数群をまとめる。反射モデルは単純な cannonball から、完全ランバート／完全鏡面、Cook–Torrance、Ashikhmin–Shirley (異方性 Phong BRDF) までを揃える。cannonball 版 (`srpCannon`) 以外は、`readSC` で読み込んだ衛星ファセット構造体 `sat` (フィールド `normal`, `area`, `pos`, `qlb`, 光学係数 `Cd`/`Cs`/`Ca`/`F0`/`kappa`/`mCT`/`nu`/`nv` 等) と、`orbitConst` の定数構造体 `const` (`S0`, `c` を参照) を受け取り、各ファセット (nFacet 行) を行方向に並べた nx3 行ベクトル規約で力・トルクを返す。共通の前段として太陽ベクトル `sunB` を正規化し、太陽距離 `d` [m] を `km2au` で AU に変換して係数 `-S0/(c·dAU^2)` を作る。数値積分版 (Cook–Torrance, Ashikhmin–Shirley) は局所座標 (法線を z 軸に一致させた系) へ `qRotation`/`qInv` で回転してからモンテカルロ積分する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `srpCannon` | cannonball モデルによる SRP 加速度 (慣性系, 面積対質量比ベース) |
| `srpSimple` | 完全ランバート反射＋完全鏡面反射モデルによる SRP 力・トルク (解析) |
| `srpCT` | Cook–Torrance モデル＋重点サンプリングによる SRP 力・トルク (数値積分) |
| `srpCTuni` | Cook–Torrance モデル＋一様サンプリングによる SRP 力・トルク (非推奨) |
| `srpCTinterp` | 補正係数テーブルの内挿で Cook–Torrance 相当の SRP を高速計算 |
| `srpAS` | Ashikhmin–Shirley モデル＋重点サンプリングによる SRP 力・トルク (数値積分) |
| `srpASuni` | Ashikhmin–Shirley モデル＋一様サンプリングによる SRP 力・トルク (非推奨) |
| `makeCoeffTable` | Cook–Torrance の補正係数テーブル `correctionPara.mat` を生成する検証用 (スクリプト) |

### Cannonball / 単純モデル

#### `srpCannon`
```matlab
aSRP = srpCannon(satAm, sunRel, d, const, Cr)
```
cannonball モデルで SRP 加速度を計算する。形状を球で近似し、面積対質量比と反射率のみで加速度を与える。
- **入力**: `satAm` (スカラ, m^2/kg) — 面積対質量比; `sunRel` (nx3, 単位ベクトル) — 衛星→太陽方向 (慣性系); `d` (nx1, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数 (`S0`, `c` を使用); `Cr` (スカラ, 無次元, optional, 既定 1.2) — 反射率係数
- **出力**: `aSRP` (nx3, m/s^2) — 慣性系での SRP 加速度
- **Note**: 内部で `sunRel` を各行正規化する。加速度出力のため質量 (`satAm`) を含む。出典 Montenbruck & Gill, *Satellite Orbits*, p79。
- *See also*: `readSC`, `orbitConst`

#### `srpSimple`
```matlab
[sat, srpCdOut, srpCsOut] = srpSimple(sat, sunB, d, const)
```
完全ランバート反射・完全鏡面反射モデルで各ファセットの SRP 力・トルクを解析的に計算する。数値積分を行わない最も軽い面ごとモデル。
- **入力**: `sat` (struct) — 衛星構成 (`normal`,`area`,`pos`,`Cd`,`Cs`,`Ca`,`kappa` 等); `sunB` (1x3, 単位ベクトル) — 衛星→太陽方向 (機体固定系); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数
- **出力**: `sat` (struct) — 入力に `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分の総和; `srpCsOut` (1x3, N) — 鏡面成分の総和
- **Note**: `%#codegen` 指定あり (コード生成対応)。`srpCsOut` は鏡面反射「力」であり単位は N (ドキュメントコメントは "Nm" と記すが実体は力の N)。
- *See also*: `readSC`, `orbitConst`

### Cook–Torrance モデル

#### `srpCT`
```matlab
[sat, srpCdOut, srpCsOut] = srpCT(sat, sunB, d, const, NDF, nMC)
```
Cook–Torrance BRDF と重点サンプリング (Beckmann 分布) でモンテカルロ積分し、SRP 力・トルクを計算する。鏡面成分を数値積分、拡散成分を解析で扱う。
- **入力**: `sat` (struct) — 衛星構成 (粗さ `mCT`, `F0` 等を使用); `sunB` (1x3, ベクトル) — 衛星→太陽方向 (機体固定系, 内部で正規化); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数; `NDF` (char, optional, 既定 `'Beckmann'`) — 法線分布関数。`'Beckmann'` 以外を渡すと Gauss 分布; `nMC` (スカラ, optional, 既定 10^4) — 積分サンプル数
- **出力**: `sat` (struct) — `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分総和; `srpCsOut` (1x3, N) — 鏡面成分総和
- **Note**: `NDF='Beckmann'` は重点サンプリングで分布関数が重みから相殺され、それ以外は一様サンプリングした反射方向に Jacobian `sin(thetaR)` が残る (`srpCTuni` と同系)。NaN 項は 0 に置換される。既定 `nMC` は実装 (`nargin<6` 分岐) では 10^4。
- *See also*: `readSC`, `orbitConst`

#### `srpCTuni`
```matlab
[sat, srpCdOut, srpCsOut] = srpCTuni(sat, sunB, d, const, NDF, nMC)
```
Cook–Torrance BRDF を全球一様サンプリングで積分する版 (非推奨)。`srpCT` と結果は等価だが収束が悪いため参照・検証用。
- **入力**: `sat` (struct) — 衛星構成; `sunB` (1x3, ベクトル) — 衛星→太陽方向 (機体固定系, 内部で正規化); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数; `NDF` (char, optional, 既定 `'Beckmann'`) — 法線分布関数。`'Beckmann'` 以外で Gauss 分布; `nMC` (スカラ, optional, 既定 10^6) — 積分サンプル数
- **出力**: `sat` (struct) — `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分総和; `srpCsOut` (1x3, N) — 鏡面成分総和
- **Note**: 半球ではなく全球 (thetaR ∈ [0,π]) を一様サンプルするため大きな `nMC` が必要 (既定 10^6)。ドキュメントコメントは第5引数を "varagin" と呼ぶが実体は名前付き引数 `NDF`。
- *See also*: `readSC`, `orbitConst`

#### `srpCTinterp`
```matlab
[sat, srpCdOut, srpCsOut] = srpCTinterp(sat, sunB, d, const, correctionPara)
```
`makeCoeffTable` で作った補正係数 (deltaS, deltaN) を粗さ mCT と入射角で内挿し、モンテカルロ積分なしに Cook–Torrance 相当の SRP を高速計算する。
- **入力**: `sat` (struct) — 衛星構成 (`mCT`, `Cd`, `Cs` 等); `sunB` (1x3, ベクトル) — 衛星→太陽方向 (機体固定系, 内部で正規化); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数; `correctionPara` (struct) — 補正係数 (`deltaS`, `deltaN`, `mSpan`, `thetaIspan`)
- **出力**: `sat` (struct) — `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分 (近似分離); `srpCsOut` (1x3, N) — 鏡面成分 (総和から拡散・入射成分を差し引いた残り)
- **Note**: `interp2` の外挿値は 1 (標準鏡面) にクランプされる。`srpCdOut`/`srpCsOut` は総力からの近似分離であり、他の関数のような厳密な成分分離ではない。ドキュメントコメントは optional の `NDF`/`nMC` を記載するが、実装のシグネチャに該当引数は無い (受け取らない)。
- *See also*: `readSC`, `orbitConst`

### Ashikhmin–Shirley モデル

#### `srpAS`
```matlab
[sat, srpCdOut, srpCsOut] = srpAS(sat, sunB, d, const, nMC)
```
Ashikhmin–Shirley の異方性 Phong BRDF と重点サンプリングで SRP 力・トルクを計算する。異方性パラメータ `nu`, `nv` に対応。拡散成分は解析、鏡面成分は数値積分。
- **入力**: `sat` (struct) — 衛星構成 (`nu`,`nv`,`F0`,`Cd`,`uu`,`uv` 等を使用); `sunB` (1x3, ベクトル) — 衛星→太陽方向 (機体固定系, 内部で正規化); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数; `nMC` (スカラ, optional, 既定 10^4) — 積分サンプル数
- **出力**: `sat` (struct) — `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分総和; `srpCsOut` (1x3, N) — 鏡面成分総和
- **Note**: 半ベクトル h を分布に沿って重点サンプリング (4象限に分割)。無限大となる重み項は 0 に置換。出典 Ashikhmin & Shirley (2000)。
- *See also*: `readSC`, `orbitConst`

#### `srpASuni`
```matlab
[sat, srpCdOut, srpCsOut] = srpASuni(sat, sunB, d, const, nMC)
```
Ashikhmin–Shirley BRDF を全球一様サンプリングで積分する版 (非推奨)。`srpAS` の参照・検証用で収束は悪い。
- **入力**: `sat` (struct) — 衛星構成 (`nu`,`nv`,`F0`,`Cd`,`uu`,`uv` 等); `sunB` (1x3, ベクトル) — 衛星→太陽方向 (機体固定系, 内部で正規化); `d` (スカラ, m) — 衛星・太陽間距離; `const` (struct) — 軌道定数; `nMC` (スカラ, optional, 既定 10^5) — 積分サンプル数
- **出力**: `sat` (struct) — `sat.force` (nx3, N), `sat.torque` (nx3, Nm) を追加; `srpCdOut` (1x3, N) — 拡散成分総和; `srpCsOut` (1x3, N) — 鏡面成分総和
- **Note**: 反射方向 v を全球 (thetaR ∈ [0,π]) 一様サンプルし Jacobian `sin(thetaR)` を重みに含める。大きな `nMC` が必要 (既定 10^5)。
- *See also*: `readSC`, `orbitConst`

### 検証用スクリプト

#### `makeCoeffTable`
```matlab
makeCoeffTable   % スクリプト（引数・戻り値なし）
```
与えた `F0` に対して、粗さ mSpan と入射角グリッドで `srpCT` (真値) と `srpSimple` を比較し、Cook–Torrance を単純モデルへ補正する係数テーブルを作る。結果を `correctionPara.mat` (deltaS, deltaN, mSpan, thetaIspan) に保存し、補正の妥当性を図示する。(スクリプト)
- **Note**: 先頭で `clear` するため実行前のワークスペース変数は消える。中心差分・2次片側差分で補正係数の勾配も計算。生成物 `correctionPara.mat` は `srpCTinterp` の入力。ユーティリティ `cls`, `gFigs`, `au2km`, `km2au`, `expint` に依存。

---

## 高精度 SRP 近似 (`hifiSRP/`)

Cook–Torrance 反射モデルを用いて、平板 facet の鏡面反射に由来する太陽輻射圧 (SRP) を解析的に近似するモジュール。可視半球を φ/θ 方向にメッシュ分割し、各セルの反射ローブを指数関数の 1 次近似で置き換えて閉形式積分することで、モンテカルロを回さずに SRP 力ベクトルを得る。衛星形状は `readSC` で読む `sat` 構造体（`sat.normal` / `sat.area` / `sat.F0` / `sat.mCT` を参照）、物理定数は `const` 構造体（`S0` / `c`）で渡す。出力の SRP 力は **太陽固定座標系** で表現される。トップの中核関数は `srpApproxCT`（単一 facet）と `srpApproxCT2`（複数 facet 対応）で、残りは係数計算・積分境界・解析解を担うヘルパー、`integration.m` は導出用スクリプト。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `srpApproxCT` | Cook–Torrance モデルによる SRP 鏡面項の解析近似（単一 facet） |
| `srpApproxCT2` | 同上の複数 facet 対応版（`sat.area` の行数だけループして総和） |
| `ctM` | Cook–Torrance モデルの NDF 以外の項 M = G·F/4 を計算（横/縦ベクトル両対応） |
| `ctM2` | 同上を facet ごと（nFacet×3）に一括計算する版 |
| `calcCoeff` | 各 φ/θ セルの 1 次近似係数 α, β を計算 |
| `intBound` | 1 セル内での鏡面ローブ内積の積分境界 [a, b] を計算 |
| `analyticSolCT` | 1 セルの解析積分 A, B（1 次近似の閉形式解）を評価 |
| `integration` | A, B の閉形式を導出したシンボリック計算ノート (スクリプト) |

### Core

#### `srpApproxCT`
```matlab
srp = srpApproxCT(sat, thetaN, sunB, d, const)
```
Cook–Torrance モデルで SRP の鏡面反射項を解析近似する（単一 facet 想定）。可視半球を +y 側の四分球と -y 側の部分半球に分けて φ/θ を 6×6 でメッシュ分割し、各セルの `calcCoeff`→`analyticSolCT` を積算する。
- **入力**: `sat` (struct) — 衛星形状（`normal` は内部で `[0, sin θN, cos θN]` に上書き、`area` / `F0` / `mCT` を使用）; `thetaN` (:,1, rad) — 太陽ベクトルと facet 法線のなす角（`arguments` は列ベクトルを要求）; `sunB` ((:,3)) — 衛星→太陽ベクトル（機体固定系）; `d` (:,1, m) — 衛星–太陽間距離; `const` (struct) — 定数（`S0` W/m^2, `c` m/s）
- **出力**: `srp` (1×3) — 近似 SRP 力ベクトル（太陽固定系）
- **Note**: 内部で `d` を `km2AU(d/1e3, const)` に通すため `d` は m 単位で渡す。φ/θ 分割数は 6×6 固定（ハードコード）。`sat.normal` は引数を無視して θN から再構成される。

#### `srpApproxCT2`
```matlab
srp = srpApproxCT2(sat, thetaN, sunB, d, const)
```
`srpApproxCT` の複数 facet 版。`sat.area` の行数 nFacet だけループし、facet ごとに `srpApproxCT` と同じメッシュ積分を行って総和した 1×3 力ベクトルを返す。
- **入力**: `sat` (struct) — `normal` は `[0, sin θN, cos θN]` の nFacet×3 に上書き、`area` / `F0` / `mCT` を使用; `thetaN` (:,1, rad) — facet ごとの太陽–法線角（nFacet 要素）; `sunB` ((:,3)) — 衛星→太陽ベクトル（機体固定系）; `d` (:,1, m) — 衛星–太陽間距離; `const` (struct) — 定数
- **出力**: `srp` (1×3) — 全 facet 合計の近似 SRP 力ベクトル（太陽固定系）
- **Note**: この版では `sunlitFlag` が `1` に固定されており、日陰判定は無効化されている（コメントアウト済み）。`srpApproxCT` は法線と太陽の内積で影を落とす点が異なる。

### Helpers

#### `ctM`
```matlab
M = ctM(sat, v, sunB)
```
Cook–Torrance モデルの NDF（法線分布関数）項を除いた残余項 M(v) = G(v)·F(v)/4 を計算する。G は幾何減衰、F は Fresnel 項。`v` と `sunB` は縦/横どちらでも受け付け内部で 3×1 に整形する。
- **入力**: `sat` (struct) — `normal` (nFacet×3), `F0` (スカラ or nFacet×1); `v` (1×3 or 3×1) — 参照ベクトル; `sunB` (1×3 or 3×1) — 衛星→太陽ベクトル（機体固定系）
- **出力**: `M` (スカラ) — 残余項（facet 方向に総和した値）
- **Note**: `n = sat.normal` を行列として扱い、複数 facet の寄与を最後に `sum(...,1)` で 1 つに畳む。2026-07-06 に `srp/ctM.m` と重複実装を統合済み。

#### `ctM2`
```matlab
M = ctM2(sat, v, sunB)
```
`ctM` を facet ごとに一括計算する版。総和せず nFacet×1 の M を返す（`srpApproxCT2` が facet ループ内で各成分を使う）。
- **入力**: `sat` (struct) — `normal` (nFacet×3), `area`, `F0` (nFacet×1); `v` (nFacet×3) — 参照ベクトル（facet ごと）; `sunB` (nFacet×3) — 衛星→太陽ベクトル（機体固定系）
- **出力**: `M` (nFacet×1) — facet ごとの残余項
- **Note**: doc コメントは `v` / `sunB` を 1×3 と記すが、実装は `vecnorm(h,2,2)` と `dot(...,2)` で **nFacet×3 行列**を前提に動く（下記 inconsistency 参照）。

#### `calcCoeff`
```matlab
[alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN, n, lam, mu)
```
メッシュ各セルについて、鏡面ローブの指数関数を線形（1 次）近似したときの係数 α, β を計算する。セル内の内積範囲 [a, b] を `intBound` で求め、その区間で exp を 1 次フィットする。
- **入力**: `phiBound` (nθ+1 × nφ+1) — φ 方向メッシュ境界（rad）; `thetaBound` (同サイズ) — θ 方向メッシュ境界（rad）; `thetaN` (スカラ, rad) — 太陽–法線角; `n` (1×3) — facet 法線ベクトル; `lam` (スカラ) — Cook–Torrance 鋭さ係数 λ = 2/mCT^2; `mu` (スカラ) — 振幅係数
- **出力**: `alp_` (nθ × nφ) — セルごとの α 係数; `bet_` (nθ × nφ) — セルごとの β 係数
- **Note**: セルが可視半球の境界（φ=π/2 かつ θN を跨ぐ）に該当するかを `nFlag` で判定し、`intBound` の `lobeIn` 引数として渡す。出力サイズは境界配列より各次元 1 小さい。

#### `intBound`
```matlab
[a, b] = intBound(phi, theta, n, lobeIn)
```
1 セル（φ×θ 矩形）4 隅の方向ベクトルと法線 `n` の内積から、鏡面ローブ内積 v·n の積分下限 a・上限 b を求める。
- **入力**: `phi` (1×2, rad) — [φ0, φ1] セル φ 範囲; `theta` (1×2, rad) — [θ0, θ1] セル θ 範囲; `n` (1×3) — facet 法線; `lobeIn` (スカラ) — ローブ内包フラグ（`1`: 上限を 1 にクランプ / `-1`: 下限を -1 にクランプし上限は 4 隅最大 / その他: 上下限とも 4 隅の min/max）
- **出力**: `a` (スカラ) — 内積の下限; `b` (スカラ) — 内積の上限
- **Note**: `lobeIn` の 3 値で境界クランプの挙動が切り替わる。ローブのピークがセルに含まれる場合の端点処理に使う。

#### `analyticSolCT`
```matlab
[Axyz, Bxyz] = analyticSolCT(thetaN, alp, bet, phi, theta)
```
1 セルの積分を、`calcCoeff` の 1 次係数 α, β を使って閉形式で評価する。`integration.m` で導出した解析式に φ/θ 境界を代入し、α 倍の項 Axyz と β 倍の項 Bxyz を返す。
- **入力**: `thetaN` (スカラ, rad) — 太陽–法線角; `alp` (スカラ) — セルの α 係数; `bet` (スカラ) — セルの β 係数; `phi` (1×2, rad) — [φ0, φ1]; `theta` (1×2, rad) — [θ0, θ1]
- **出力**: `Axyz` (1×3) — α 項の積分寄与（α でスケール済み）; `Bxyz` (1×3) — β 項の積分寄与（β でスケール済み）
- **Note**: 式は半角（θ/2）を多用した長大な多項式で、`integration.m` のシンボリック導出が出典。呼び出し側でセルごとに積算する。

### Scripts

#### `integration`
```matlab
% スクリプト（引数・戻り値なし）
```
`analyticSolCT` の閉形式解を導いた Symbolic Math Toolbox のライブスクリプト。半球ベクトル h、法線 n、Cook–Torrance の指数重み σ1 を定義し、被積分関数を φ/θ で不定積分して A ベクトルを得る導出過程を記録する。
- **Note**: 冒頭で `clc; cls; clear` を実行する開発用ノートで、関数としては呼べない（`cls` は標準 MATLAB コマンドではなく環境依存）。実運用コードではなく導出の再現用。

---

## ライトカーブ (`lightcurves/`)

衛星の光度曲線（ライトカーブ）を計算する関数群。`object/readSC` で読み込んだ衛星モデル構造体 `sat`（`normal` / `area` / `Cd` / `Cs` / `F0` / `nu` / `nv` / `uu` / `uv` / `mCT` 等の面ごとの反射パラメータを持つ）を入力とし、指定した BRDF（双方向反射率分布関数）モデルで各面の反射光量を積算する。時系列/観測方向は横ベクトル規約で扱われ、多くの中間量は `nFacet × M`（面数 × 時刻数）の行列になる。メイン関数 `lc` が観測ジオメトリ（慣性系の位置ベクトルと姿勢クォータニオン）から body-fixed 系の太陽・観測者方向を作り、`lcSimple` / `lcAS` / `lcCT` のいずれかを呼び出して見かけの等級を返す。反射光量と等級の相互変換は `mag` / `magInv`、モデルの反射特性の可視化は `visualizeBRDF`（スクリプト）。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `lc` | ライトカーブ計算のメイン関数。姿勢・位置からBRDFモデルを選んで見かけ等級を返す |
| `lcSimple` | Lambertian拡散＋完全鏡面反射の簡易BRDFモデル |
| `lcAS` | Ashikhmin–Shirleyモデル（異方性反射）による反射光量計算 |
| `lcCT` | Cook–Torranceモデル（面粗さ考慮）による反射光量計算 |
| `mag` | 反射光量と距離から見かけの等級を計算 |
| `magInv` | 見かけの等級と距離から反射光量を逆算 |
| `visualizeBRDF` | BRDFの反射特性を半球面上に3Dプロットで可視化 (スクリプト) |

### Core

#### `lc`
```matlab
[m, fObs] = lc(sat, scalar, q, satPos, obsPos, sunPos, nu, options)
```
ライトカーブ計算のメイン関数。慣性系の位置ベクトルと姿勢から body-fixed 系での太陽・観測者方向を作り、`options.BRDF` で選んだモデルで反射光量を積算して太陽基準の相対等級を返す。
- **入力**: `sat` (struct) — 衛星モデル; `scalar` (1×1) — クォータニオンのスカラー部位置（`0`=先頭 [cos(θ/2), e·sin(θ/2)]、`4`=末尾 [e·sin(θ/2), cos(θ/2)]）; `q` (n×4) — 姿勢クォータニオン; `satPos` (n×3, m) — 慣性系の衛星位置; `obsPos` (n×3, m) — 慣性系の観測者位置; `sunPos` (n×3, m) — 慣性系の太陽位置; `nu` (n×1) — 日照フラグ（`1`=sunlit, `0`=eclipse）; `options` (struct, 任意) — `options.BRDF`（`'simple'`(既定) / `'AS'` / `'CT'`）と `options.mex`（`'on'`/`'off'`(既定)）
- **出力**: `m` (n×1) — 太陽基準の見かけ相対等級の時系列; `fObs` (n×1) — 全面積算した反射光量の時系列（日照フラグ適用済み）
- **Note**: BRDF モデルは位置引数ではなく `options.BRDF` フィールドで指定する。等級は `m = -26.7 - 2.5·log10(fObs / |obsPos-satPos|^2)` で計算され、距離は m 単位。`options.mex = 'on'` の場合は `lcSimple_mex` 等の生成コードを呼ぶ。
- *See also*: `readSC`, `lcSimple`, `lcAS`, `lcCT`, `orbitConst`

### Variants（BRDFモデル）

これらは `lc` から呼ばれるほか単独でも使える。入力の `sunB` / `obsB` は body-fixed 系での単位方向ベクトルで、返り値の主産物は `sat.fObs`（`nFacet × M`）フィールドに格納される。

#### `lcSimple`
```matlab
[sat, cd, cs] = lcSimple(sat, sunB, obsB, thr)
```
Lambertian 拡散反射と完全鏡面反射を組み合わせた簡易 BRDF モデル。鏡面成分は太陽・観測者の二等分線が面法線と閾値 `thr` 内で一致したとき点灯する。
- **入力**: `sat` (struct) — 衛星モデル; `sunB` (M×3) — body-fixed 系での衛星→太陽の単位ベクトル; `obsB` (M×3) — body-fixed 系での衛星→観測者の単位ベクトル; `thr` (1×1, rad, 任意) — 鏡面ローブの閾値（既定 `deg2rad(1)`）
- **出力**: `sat` (struct) — `sat.fObs`（nFacet×M）が追加された衛星モデル; `cd` (nFacet×M) — 拡散反射成分; `cs` (nFacet×M) — 鏡面反射成分
- **Note**: `sunB` / `obsB` は内部で正規化される。面が可視（`NS>0` かつ `NV>0`）でない要素は 0 にマスクされる。
- *See also*: `lcAS`, `lcCT`, `readSC`

#### `lcAS`
```matlab
[sat, cd, cs, D] = lcAS(sat, sunB, obsB)
```
Ashikhmin–Shirley 異方性 Phong BRDF モデルによる反射光量計算。面上の直交接ベクトル `sat.uu` / `sat.uv` と指数 `sat.nu` / `sat.nv` で異方性を表現する。
- **入力**: `sat` (struct) — 衛星モデル（N facets、`normal` / `area` / `Cd` / `F0` / `nu` / `nv` / `uu` / `uv` を使用）; `sunB` (M×3) — body-fixed 系での衛星→太陽の単位ベクトル; `obsB` (M×3) — body-fixed 系での衛星→観測者の単位ベクトル
- **出力**: `sat` (struct) — `sat.fObs`（nFacet×M）が追加された衛星モデル; `cd` (nFacet×M) — 拡散反射成分; `cs` (nFacet×M) — 鏡面反射成分; `D` (nFacet×M) — 異方性分布項
- **Note**: body-fixed 系で計算するため面法線は必ずしも [0,0,1] ではない。`sunB` / `obsB` は内部で正規化される。`max(NS,NV)` による 0 割で生じる Inf は 0 に置換され、複素数化を防ぐため `NS<0` / `NV<0` の要素は先にマスクされる。
- *See also*: `readSC`, `lcSimple`, `lcCT`

#### `lcCT`
```matlab
[sat, cd, cs, D] = lcCT(sat, sunB, obsB, NDF)
```
Cook–Torrance モデルによる物理ベースの反射光量計算。微小面の粗さ `sat.mCT`、幾何減衰項 G、フレネル項 F を考慮する。
- **入力**: `sat` (struct) — 衛星モデル（`normal` / `area` / `Cd` / `F0` / `mCT` を使用）; `sunB` (M×3) — body-fixed 系での衛星→太陽の単位ベクトル; `obsB` (M×3) — body-fixed 系での衛星→観測者の単位ベクトル; `NDF` (char, 任意) — 法線分布関数 `'Beckmann'`(既定) または `'Gauss'`
- **出力**: `sat` (struct) — `sat.fObs`（nFacet×M）が追加された衛星モデル; `cd` (nFacet×M) — 拡散反射成分; `cs` (nFacet×M) — 鏡面反射成分; `D` (nFacet×M) — 法線分布関数
- **Note**: `NDF` に `'Gauss'` 以外の文字列を渡すと `'Beckmann'` にフォールバックする。`sunB` / `obsB` は内部で正規化される。面が可視でない要素は 0 にマスクされる。
- *See also*: `lcAS`, `readSC`

### Helpers

#### `mag`
```matlab
appMag = mag(fobs, d)
```
反射光量と距離から見かけの等級を計算する。太陽の見かけ等級 -26.7 を基準にする。
- **入力**: `fobs` — 反射光量（フラックス）; `d` (m) — 観測者と対象の距離（`fobs` と要素ごとに演算できる形状）
- **出力**: `appMag` — 見かけの等級
- **Note**: `appMag = -26.7 - 2.5·log10(fobs / d^2)`。距離は m 単位。

#### `magInv`
```matlab
fobs = magInv(mApp, d)
```
見かけの等級と距離から反射光量を逆算する（`mag` の逆演算）。
- **入力**: `mApp` (N×1) — 太陽基準の見かけ相対等級; `d` (N×1, m) — 観測者と対象の距離
- **出力**: `fobs` (N×1) — 反射光量（フラックス）
- **Note**: `fobs = d^2 · 10^((mApp+26.7)/(-2.5))`。距離は m 単位。
- *See also*: `readSC`, `lcSimple`, `lcAS`, `lcCT`, `orbitConst`

### Scripts

#### `visualizeBRDF`
```matlab
visualizeBRDF   % スクリプト（引数なし）
```
単一面の `sat` パラメータと固定した太陽入射方向に対し、半球面上の反射強度分布を `surf` で 3D 可視化するスクリプト。`lcAS`（既定）または `lcCT`（コメントアウト）を反射方向のグリッド上で評価する。
- **Note**: 先頭で `clc` / `clear` / `cls` を実行するためワークスペースをクリアする。`sat.nu` / `sat.nv` / `sat.uu` / `sat.uv`（AS 用）、`sat.mCT`（CT 用）等をスクリプト内で直接設定する。

---

## 球面ガウス関数 (`sphericalGaussian/`)

球面ガウス分布 (Spherical Gaussian, SG) とその混合の計算を扱う小さな関数群。SG は単位球面上の点 x に対し μ·exp(λ(x·p − 1)) で定義される lobe 型の分布で、反射・照明モデル等の基底として用いる。座標は横ベクトル規約ではなく、球面上の座標を x/y/z 各成分の配列（`meshgrid` 由来の NxM 行列でも可）として個別に渡す設計になっている。lobe 方向 p は 1x3 列/行ベクトル、sharpness λ・amplitude μ はスカラー（`sgMix` のみ n 個をまとめた縦ベクトルも可）。外部構造体（`sat`/`const`/`oe` 等）には依存しない。可視化用のデモスクリプト `visualizeSG` もここに含む。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `sg` | 等方球面ガウス μ·exp(λ(x·p−1)) を球面座標上で評価する |
| `asg` | 異方性球面ガウス（tangent/bi-tangent/lobe 軸で帯域幅を分ける）を評価する |
| `sgMix` | 2 つの SG の積を 1 つの SG（p, λ, μ）に閉じた形で合成する |
| `visualizeSG` | SG を球面上に surf 描画するデモ (スクリプト) |

#### `sg`
```matlab
G = sg(x, y, z, p, lam, mu)
```
等方な球面ガウス μ·exp(λ(x·p−1)) を、球面上の各座標について評価する。関数内で lobe 軸 p は単位ベクトルに正規化される。
- **入力**: `x` (NxM) — 球面上の x 座標; `y` (NxM) — y 座標; `z` (NxM) — z 座標; `p` (1x3) — lobe 方向ベクトル（内部で正規化）; `lam` (スカラー) — sharpness λ; `mu` (スカラー) — amplitude μ
- **出力**: `G` (x,y,z と同形) — SG の値。doc は「Nx1 vector」と書くが、実装は入力配列と同じ形状の配列を返す
- *See also*: （doc の `sgMult`, `sg3` はライブラリに存在しない）

#### `asg`
```matlab
aSG = asg(vx, vy, vz, x, y, z, lam, mu, c)
```
異方性球面ガウス。球面上の点を tangent 軸 x・bi-tangent 軸 y・lobe 軸 z へ射影し、lobe 軸方向の可視項 max(0, VZ) に、x/y 各方向で帯域幅の異なるガウス減衰 exp(−lam·VX² − mu·VY²) を掛ける。
- **入力**: `vx, vy, vz` (各 NxM) — 球面上の点の x/y/z 成分; `x` (1x3) — tangent 軸; `y` (1x3) — bi-tangent 軸; `z` (1x3) — lobe 軸; `lam` (スカラー) — x 軸方向の帯域幅; `mu` (スカラー) — y 軸方向の帯域幅; `c` (スカラー) — amplitude
- **出力**: `aSG` (入力配列と同形) — 異方性 SG の値
- **Note**: `mu` は amplitude ではなく y 軸方向の帯域幅であり、amplitude は別引数 `c`（等方 `sg` の μ とは役割が異なる点に注意）。lobe の裏側は max(0, VZ) により 0 にクランプされる。
- *See also*: （doc の `sgMult` はライブラリに存在しない）

#### `sgMix`
```matlab
[p3, lam3, mu3] = sgMix(p1, lam1, mu1, p2, lam2, mu2)
```
2 つの球面ガウスの積を 1 つの SG に閉じた形で合成する。lam3 = ‖λ1·p1 + λ2·p2‖、p3 = (λ1·p1 + λ2·p2)/lam3、mu3 = μ1·μ2·exp(lam3 − (λ1 + λ2))。
- **入力**: `p1` (nx3) — SG1 の lobe 軸; `lam1` (nx1) — SG1 の sharpness; `mu1` (nx1) — SG1 の amplitude; `p2` (nx3) — SG2 の lobe 軸; `lam2` (nx1) — SG2 の sharpness; `mu2` (nx1) — SG2 の amplitude
- **出力**: `p3` (nx3) — 合成 SG の lobe 軸; `lam3` (nx1) — 合成 sharpness; `mu3` (nx1) — 合成 amplitude
- **Note**: `arguments` ブロックは各引数を `(:,1)`/`(:,3)` の縦ベクトル/n行として宣言しており、複数 SG をまとめて一括合成できる（doc の「scalar」表記より一般）。
- *See also*: `sg`

#### `visualizeSG` (スクリプト)
SG の挙動を確認するためのデモスクリプト。`meshgrid` で球面グリッドを作り、`sg` で SG 値を求め、`surf` で球面上に色付き描画し lobe 軸を `quiver3` で重ねる。引数・戻り値はなく、冒頭で `clc/clear` する。パラメータ（lobe 軸 p1/p2、λ、μ）はスクリプト内にハードコードされている。

---

## UKF / CKF フィルタ (`ukfCkf/`)

Unscented Kalman Filter (UKF)、Cubature Kalman Filter (CKF)、および球面放射 Square-Root CKF (SR-CKF) の構成要素関数群。各フィルタは「シグマ点生成 → 予測共分散 → 相互共分散とカルマンゲイン」の3ステップに分かれ、それぞれ `*Sigma` / `*Cov` / `*CorrGain` 関数が対応する。ベクトル規約は横ベクトル（行方向）が主で、シグマ点行列は「行=各シグマ点、列=状態次元」の `(点数) x n` 形式（UKF/SR-CKF は 2n+1 点、CKF は 2n 点）。関数群は素の行列引数で動き、専用の構造体依存はない（`ukfInitPara` のみ `ukf_` 構造体にチューニングパラメータ・重みを格納する）。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `ukf` | UKF の1ステップ（予測・更新）をまとめて実行する統合関数 |
| `ukfSigma` | 共分散平方根（SVD/Cholesky）からシグマ点 (2n+1)×n を生成 |
| `ukfCov` | シグマ点から事前共分散（プロセスノイズ Q 加算）を計算 |
| `ukfCorrGain` | 観測共分散 Pyy・相互共分散 Pxy・カルマンゲイン K を計算 |
| `ukfInitPara` | UKF のチューニングパラメータ（α, β, κ, λ）と重み wm/wc を初期化 |
| `ckfSigma` | Cholesky 平方根から球面キュバチャ点 2n×n を生成（中心点なし） |
| `ckfCov` | キュバチャ点から事前共分散を計算（等重み 1/2n） |
| `ckfCorrGain` | CKF の Pyy・Pxy・K を計算（等重み） |
| `srckfSigma` | 中心点付き球面放射キュバチャ点 (2n+1)×n を生成 |
| `srckfCov` | 重み wc を用いた SR-CKF 事前共分散を計算 |
| `srckfCorrGain` | 重み wc を用いた SR-CKF の Pyy・Pxy・K を計算 |

### Core (UKF)

#### `ukf`
```matlab
[Pest, xEst] = ukf(f, h, Q, R, lam, P, x, y, wm, wc, wic)
```
UKF の1サイクル（シグマ点生成・伝播・予測・観測更新）を実行する統合関数。
- **入力**: `f` (関数ハンドル) — 離散システムダイナミクス x(k+1)=f(x); `h` (関数ハンドル) — 観測モデル y=h(x); `Q` (nxn) — プロセスノイズ共分散; `R` (mxm) — 観測ノイズ共分散; `lam` (スカラー) — UKF チューニングパラメータ λ; `P` (nxn) — 事前共分散; `x` (nx1) — 状態推定値; `y` (mx1) — 観測ベクトル; `wm`/`wc`/`wic` — 平均・共分散用の重み
- **出力**: `Pest` (nxn) — 更新後共分散; `xEst` (1xn) — 更新後状態推定値
- **Note**: 実装に不具合あり。31 行目 `xEst = sum(wim .* X, 1)` が引数リストに存在しない変数 `wim` を参照しており、そのままでは実行時エラーになる（`wm` の誤記と思われる）。使用前に要修正。

#### `ukfSigma`
```matlab
X = ukfSigma(lam, P, x, method)
```
共分散行列の平方根からアンセンテッド変換のシグマ点を生成する。
- **入力**: `lam` (スカラー) — UKF チューニングパラメータ λ; `P` (nxn) — 共分散行列; `x` (nx1 または 1xn) — 状態ベクトル; `method` ('svd' または 'chol'、既定 'svd') — 平方根計算法（SVD または Cholesky 分解）
- **出力**: `X` ((2n+1)×n) — シグマ点（1 行目が中心点 x、以降が ±sqrt(n+λ)·Psq のオフセット）
- **Note**: `method` 未指定時は SVD を使用。'svd'/'chol' 以外を渡すとエラー。
- *See also*: `ukfCov`

#### `ukfCov`
```matlab
Pcov = ukfCov(xEst, X, wc, Q)
```
シグマ点と共分散重みから事前共分散を計算し、プロセスノイズ Q を加算する。
- **入力**: `xEst` (1xn) — 事前状態推定値; `X` ((2n+1)×n) — シグマ点; `wc` (スカラーまたは (2n+1) ベクトル) — 共分散用重み; `Q` (nxn) — プロセスノイズ共分散
- **出力**: `Pcov` (nxn) — 事前共分散
- *See also*: `ukfSigma`, `ukf`

#### `ukfCorrGain`
```matlab
[Pyy, Pxy, K] = ukfCorrGain(xEst, X, yEst, Y, wc, R)
```
観測共分散 Pyy、状態-観測相互共分散 Pxy、カルマンゲイン K を計算する。
- **入力**: `xEst` (nx1) — 事前状態推定値; `X` ((2n+1)×n) — 状態シグマ点; `yEst` (mx1) — 事前観測推定値; `Y` ((2n+1)×m) — 観測シグマ点; `wc` (スカラーまたは (2n+1) ベクトル) — 共分散用重み; `R` (mxm) — 観測ノイズ共分散
- **出力**: `Pyy` (mxm) — 観測共分散; `Pxy` (nxm) — 相互共分散; `K` (nxm) — カルマンゲイン K = Pxy·Pyy^-1

#### `ukfInitPara`
```matlab
ukf_ = ukfInitPara(n_, ukf_)
```
UKF のチューニングパラメータと重みを初期化して `ukf_` 構造体に格納する。各フィールドが未設定の場合のみ既定値（α=1e-4, β=2, κ=3-n, λ=α²(n+κ)-n）を計算し、平均・共分散重み wm/wc を導出する。
- **入力**: `n_` (スカラー) — 状態次元; `ukf_` (構造体) — 既存設定を保持する構造体（該当フィールドがあればそれを尊重）
- **出力**: `ukf_` (構造体) — `n_`, `alp`, `bet`, `kappa`, `lambda`, `wm` ((2n+1)×1), `wc` ((2n+1)×1) を設定した構造体
- **Note**: 既に設定済みのフィールドは上書きしない（`isfield` で分岐）。重みは列ベクトルに正規化される。

### Variants (CKF)

#### `ckfSigma`
```matlab
Xout = ckfSigma(P, x)
```
Cholesky 平方根から球面キュバチャ点（中心点なしの 2n 点）を生成する。
- **入力**: `P` (nxn) — 共分散行列; `x` (nx1 または 1xn) — 状態ベクトル
- **出力**: `Xout` (2n×n) — キュバチャ点（±sqrt(n)·S を x に加えたもの）
- **Note**: Cholesky 分解を用いるため P は正定値であること。
- *See also*: `ckfCov`

#### `ckfCov`
```matlab
Pout = ckfCov(xEst, X, Q)
```
キュバチャ点から等重み（1/2n）で事前共分散を計算し、プロセスノイズ Q を加算する。
- **入力**: `xEst` (1xn) — 事前状態推定値; `X` (2n×n) — キュバチャ点; `Q` (nxn) — プロセスノイズ共分散
- **出力**: `Pout` (nxn) — 事前共分散
- **Note**: `Pout = mean(X'X) - xEst'·xEst + Q` の形で計算するため `xEst` は行ベクトル (1xn) 前提。
- *See also*: `ukfSigma`, `ukf`

#### `ckfCorrGain`
```matlab
[Pyy, Pxy, K] = ckfCorrGain(xEst, X, yEst, Y, R)
```
CKF の観測共分散 Pyy、相互共分散 Pxy、カルマンゲイン K を等重みで計算する。
- **入力**: `xEst` (nx1) — 事前状態推定値; `X` (2n×n) — 状態キュバチャ点; `yEst` (mx1) — 事前観測推定値; `Y` (2n×m) — 観測キュバチャ点; `R` (mxm) — 観測ノイズ共分散
- **出力**: `Pyy` (mxm) — 観測共分散; `Pxy` (nxm) — 相互共分散; `K` (nxm) — カルマンゲイン
- **Note**: UKF 版と異なり重み引数を取らず、常に 1/2n の等重みを用いる。

### Variants (Square-Root CKF)

#### `srckfSigma`
```matlab
Xout = srckfSigma(P, x)
```
中心点を含む球面放射キュバチャ点 (2n+1 点) を Cholesky 平方根から生成する。
- **入力**: `P` (nxn) — 共分散行列; `x` (nx1 または 1xn) — 状態ベクトル
- **出力**: `Xout` ((2n+1)×n) — キュバチャ点（1 行目が中心点 x、以降が x±sqrt(n)·S）
- **Note**: 中心点を含むため点数は 2n+1（中心点なしの `ckfSigma` は 2n）。ドキュメントコメントの出力次元「2n×n」は誤りで、実際は (2n+1)×n。
- *See also*: `ckfCov`

#### `srckfCov`
```matlab
Pout = srckfCov(xEst, X, wc, Q)
```
重み wc を用いて SR-CKF の事前共分散を計算し、プロセスノイズ Q を加算する。
- **入力**: `xEst` (nx1 または 1xn) — 事前状態推定値; `X` ((2n+1)×n) — キュバチャ点; `wc` ((2n+1) ベクトル) — 共分散用重み; `Q` (nxn) — プロセスノイズ共分散
- **出力**: `Pout` (nxn) — 事前共分散
- **Note**: 点数 `nSigma = size(X,1)` を用いてループするため実際の点数に追従する。ドキュメントの `X: 2n×n` は (2n+1)×n の誤り。
- *See also*: `ukfSigma`, `ukf`

#### `srckfCorrGain`
```matlab
[Pyy, Pxy, K] = srckfCorrGain(xEst, X, yEst, Y, wc, R)
```
重み wc を用いて SR-CKF の観測共分散 Pyy、相互共分散 Pxy、カルマンゲイン K を計算する。
- **入力**: `xEst` (nx1) — 事前状態推定値; `X` ((2n+1)×n) — 状態キュバチャ点; `yEst` (mx1) — 事前観測推定値; `Y` ((2n+1)×m) — 観測キュバチャ点; `wc` ((2n+1)×1) — 共分散用重み; `R` (mxm) — 観測ノイズ共分散
- **出力**: `Pyy` (mxm) — 観測共分散; `Pxy` (nxm) — 相互共分散; `K` (nxm) — カルマンゲイン K = Pxy·Pyy^-1
- **Note**: ドキュメントの `X: 2n×n`, `Y: 2n×m` は (2n+1)×n / (2n+1)×m の誤り（`wc` の記載は 2n+1×1 で正しい）。

---

## ガウス過程回帰 (GPR) (`gpr/`)

ガウス過程回帰 (GPR) の予測平均・予測分散と、その基礎となるガウスカーネルを計算する関数群。学習データ・テストデータはいずれも行方向にサンプル (データ点)、列方向に入力次元を並べる横ベクトル規約 (`xTrain` は d×n、`xAst` は m×n)。カーネル行列の逆行列を陽に扱う実装 (`gprCov`) と、数値的に安定な Cholesky 分解 (下三角 `L`) を用いる実装 (`gprMean`) が混在する点に注意。外部構造体 (`sat` / `const` / `oe` 等) への依存はない。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `gprMean` | GPR の予測平均と対数周辺尤度を計算 (Cholesky 因子 L を使用) |
| `gprCov` | GPR の予測分散をテストデータごとに計算 (カーネル逆行列を使用) |
| `kernelGauss` | 2 つの入力ベクトル間のガウスカーネル値を計算 |
| `kernelGaussMat` | 学習データ全体のガウスカーネル行列 (ノイズ項込み) を計算 |

### Core

#### `gprMean`
```matlab
[yPred, logP] = gprMean(xAst, xTrain, xMean, yTrain, yMean, L, hypPara)
```
GPR の予測平均を計算する。カーネル逆行列を陽に持つ代わりに Cholesky 因子 L を用いて α = L' \ (L \ (yTrain − xMean)) を解き、数値的に安定・高速に平均と対数周辺尤度を求める。
- **入力**: `xAst` (m×n) — テストデータ; `xTrain` (d×n) — 学習入力データ; `xMean` (d×出力次元 と整合) — 学習入力に対する平均関数の値 (yTrain から差し引く事前平均); `yTrain` (d×出力次元) — 学習出力データ; `yMean` — テスト点での事前平均 (現実装では yPred の基準として `xMean` を使用しており本文中で直接は参照されない); `L` (d×d, 下三角) — カーネル行列の Cholesky 因子; `hypPara` (1×2) — ガウスカーネルのハイパーパラメータ [θ1, θ2]
- **出力**: `yPred` (m×出力次元) — テストデータに対する予測平均; `logP` (1×1) — 対数周辺尤度 log p(y|X)
- **Note**: ドキュメントコメントの入力欄が実装と大きく食い違う (下記 inconsistencies 参照)。実装は 7 引数だが doc は 5 個 (`Kinv` 等) しか列挙していない。yPred の基準値には `yMean` ではなく `xMean` が加算されている点にも注意。
- *See also*: `gprCov`, `kernelGauss`, `kernelGaussMat`

#### `gprCov`
```matlab
sigPred = gprCov(xAst, xTrain, Kinv, hypPara)
```
GPR の予測分散 (事後共分散の対角成分) を各テストデータ点について計算する。k** − diag(k*' · Kinv · k*) を返す。
- **入力**: `xAst` (m×n) — テストデータ; `xTrain` (d×n) — 学習入力データ; `Kinv` (d×d) — カーネル行列の逆行列; `hypPara` (1×2) — ガウスカーネルのハイパーパラメータ [θ1, θ2]
- **出力**: `sigPred` (m×1) — 各テストデータ点での予測分散
- **Note**: 出力は予測平均ではなく予測分散。ドキュメントコメントの出力欄が `yPred`「予測平均」となっているのは誤り (下記 inconsistencies 参照)。doc の入力欄には実装に存在しない `yTrain` が含まれる。
- *See also*: `gprMean`, `kernelGauss`

### Helpers

#### `kernelGauss`
```matlab
k = kernelGauss(x, xp, hypPara)
```
2 つの入力ベクトル x, x' 間のガウス (RBF) カーネル値を計算する。k(x, x'; θ1, θ2) = θ1 · exp( −‖x − x'‖² / θ2 )。
- **入力**: `x` (1×n) — 入力ベクトル; `xp` (1×n) — 入力ベクトル x'; `hypPara` (1×2) — ハイパーパラメータ [θ1 (スケール), θ2 (長さスケール)]
- **出力**: `k` (1×1) — スカラーのカーネル値
- **Note**: ドキュメントコメントの出力欄が `yPred`「予測平均, Mxn matrix」となっているが、実際の出力はスカラーのカーネル値 (他関数からの copy-paste 誤り、下記 inconsistencies 参照)。
- *See also*: `gprCov`, `kernelGaussMat`

#### `kernelGaussMat`
```matlab
K = kernelGaussMat(xTrain, hypPara, sigN)
```
学習データ全体のガウスカーネル行列を計算し、対角に観測ノイズ分散を加える。K = θ1 · exp(−(ペア間距離)² / θ2) + sigN² · I。`squareform(pdist(...))` を使うため Statistics and Machine Learning Toolbox が必要 (未使用時用のループ実装がコメントで併記)。
- **入力**: `xTrain` (d×n) — 学習入力データ (行方向にデータ点); `hypPara` (1×2) — ハイパーパラメータ [θ1, θ2]; `sigN` (1×1) — 加法性ノイズの標準偏差
- **出力**: `K` (d×d) — ノイズ項込みのガウスカーネル行列
- **Note**: `squareform`/`pdist` 依存で Statistics and Machine Learning Toolbox が必須。`gather` を使っており gpuArray 入力にも対応。
- *See also*: `gprCov`, `kernelGauss`

---

## 汎用ユーティリティ (`utility/`)

この節は図の体裁調整・保存、可視化補助、簡易なプログラミングヘルパーなど、特定の力学モデルに依存しない汎用関数をまとめる。時系列データを扱う関数（`plotStd`, `drawShadowZones`）は行方向 = 時系列サンプル（N×1 の縦ベクトル）を前提とする。`sat` / `const` / `oe` 等の構造体には依存しない。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `cls` | 開いている figure ウィンドウを全て閉じる |
| `fig4Paper` | 1 枚の figure を論文・会議録用に整形して PDF 保存 |
| `figs4Paper` | 表示中の全 figure を論文用に整形して一括保存 |
| `fig4Presen` | figure をプレゼン用に整形して PDF 保存 |
| `gFigs` | 開いている figure を指定ディスプレイ上にグリッド配置 |
| `plotStd` | 平均値に標準偏差の帯を重ねてプロット |
| `drawShadowZones` | sunlit/shadow フラグに応じて eclipse 区間を xregion で網掛け |
| `sb` | symbolic 数式を Scrapbox 用 LaTeX 文字列に変換 |
| `ifelse` | 条件式による三項演算子相当のインライン分岐 |

### 図の整形・保存

#### `fig4Paper`
```matlab
fig4Paper(asis, fig, nFig, contentType)
```
論文・会議録用に figure のフォント（Times New Roman）、線幅、グリッド、LaTeX インタープリタを整えて PDF (600 dpi) 保存する。ファイル名は `fig<nFig>_<タイムスタンプ>.pdf`。ヘルパー関数（`adjustFont` / `optimizeFig` / `applyLatexLabels` / `detectLayout` / `decideFigSize`）を同一ファイル内に持つ。
- **入力**: `asis` (1×1) — 1 で figure サイズを保持、0 でレイアウトから自動リサイズ（既定 1）; `fig` — 対象 figure ハンドル（既定 `gcf`）; `nFig` (1×1) — 出力ファイル名に使う図番号（既定 1）; `contentType` — `exportgraphics` の ContentType、`"vector"` / `"image"` 等（既定 `"vector"`）
- **出力**: なし（PDF ファイルを書き出す）
- **Note**: `asis == 0` のときのみ subplot / tiledlayout の行列数を検出して図サイズ（cm 単位、1〜2 カラム幅）を決める。全 4 引数がオプションで `nargin` により既定値を補う。
- *See also*: `fig4Presen`, `figs4Paper`

#### `figs4Paper`
```matlab
figs4Paper(asis, contentType)
```
表示中（Visible）の全 figure を図番号昇順に走査し、各図に `fig4Paper` を適用して一括保存する。Live Editor のインライン figure（Tag が `TMWLiveEditor`）は除外し、`UserData` に保存済みフラグを記録して重複保存を防ぐ。
- **入力**: `asis` (1×1) — `fig4Paper` に渡すサイズ保持フラグ（既定 0）; `contentType` — `fig4Paper` に渡す ContentType（既定 `'vector'`）
- **出力**: なし（各 figure を PDF 保存）
- **Note**: `nargin` による既定値設定に不備があり、引数 0 個で呼ぶと `contentType` が未定義のまま `fig4Paper` に渡り実行時エラーになる。安全に使うには `contentType` まで含めて明示指定するか、下記の inconsistency を参照して修正すること。
- *See also*: `fig4Paper`, `fig4Presen`

#### `fig4Presen`
```matlab
fig4Presen
```
現在の座標軸 (`gca`) をプレゼン用にフォント 24pt・Times、線幅 1.5、7 色の ColorOrder、XYZ グリッド ON に設定し、`fig1.pdf`（ContentType `image`）として保存する。
- **入力**: なし（`gca` / `gcf` を操作）
- **出力**: なし（`fig1.pdf` を書き出す）
- **Note**: 出力ファイル名が `fig1.pdf` 固定で、連続実行すると上書きされる。
- *See also*: `fig4Paper`

#### `cls`
```matlab
cls
```
`close all` を実行して開いている全 figure ウィンドウを閉じるだけの短縮コマンド。
- **入力**: なし
- **出力**: なし
- *See also*: `fig4Presen`

#### `gFigs`
```matlab
gFigs(nDisp)
```
表示中の figure を指定ディスプレイ上に自動計算したグリッド（ほぼ正方形）で並べ替える。`nDisp == 0` のときは何もしない。
- **入力**: `nDisp` (1×1, 正整数) — 配置先ディスプレイ番号（1: プライマリ, 2: セカンダリ …、既定 1）
- **出力**: なし（各 figure の Position を更新）
- **Note**: `MonitorPositions`（ピクセル単位）に基づき配置するため、複数モニタ環境を前提とする。`arguments` ブロックは実装ではコメントアウトされており、既定値 1 は実際には効かず、引数なし呼び出しはエラーになる。

### 可視化補助

#### `plotStd`
```matlab
plotStd(t_, meanVal, stdVal, lineColor, regionColor)
```
平均値の折れ線に標準偏差 ±1σ の帯（半透明 area）を重ねてプロットする。`area` の第 2 引数 `[meanVal - stdVal, 2*stdVal]` を積み上げることで下端 `meanVal - stdVal`、上端 `meanVal + stdVal` の帯を描く。
- **入力**: `t_` (N×1, numeric) — x 軸の値（時刻等）; `meanVal` (N×1, numeric) — 平均値; `stdVal` (N×1, numeric) — 標準偏差; `lineColor` (text) — 折れ線の色（既定 `'r'`）; `regionColor` (text) — 帯の色（既定 `'k'`）
- **出力**: なし（現在の座標軸にプロット）
- **Note**: ドキュメントコメントは 3 入力しか記載していないが、実装は `lineColor` / `regionColor` を含む 5 入力（下記 inconsistency 参照）。`hold on` を内部で呼ぶ。
- *See also*: `fig4Paper`

#### `drawShadowZones`
```matlab
drawShadowZones(t_, sunlitFlag)
```
sunlit/shadow フラグ列に応じて、現在の図に eclipse（shadow）区間を `xregion` で網掛けする。`sunlitFlag` を両端 1 でパディングして差分を取り、1→0 を eclipse 開始、0→1 を終了として区間を検出する。全て sunlit なら何もしない。
- **入力**: `t_` (N×1) — 時刻ベクトル; `sunlitFlag` (N×1) — 1: sunlit, 0: shadow のフラグ列
- **出力**: なし（現在の座標軸に xregion を追加）
- **Note**: ドキュメント用コメントブロックを持たず、関数先頭のインラインコメントのみ。`xregion` は MATLAB R2023a 以降の組み込み関数。

### プログラミングヘルパー

#### `sb`
```matlab
out = sb(symEq)
```
symbolic 数式を `latex()` で変換し、Scrapbox 記法の `[$ ... ]` で囲んだ文字列を生成してコマンドウィンドウに出力する。
- **入力**: `symEq` — symbolic 式
- **出力**: `out` (1×1 cell) — `'[$ <latex> ]'` 形式の文字列を格納した cell
- **Note**: 出力を cell (`out{:}`) として組み立てる。表示は `fprintf` で行う。

#### `ifelse`
```matlab
result = ifelse(condition, trueValue, falseValue)
```
条件が真なら `trueValue`、偽なら `falseValue` を返すインライン三項分岐ヘルパー。
- **入力**: `condition` — 論理条件; `trueValue` — 真のときの値; `falseValue` — 偽のときの値
- **出力**: `result` — 選択された値
- **Note**: 両分岐の値は呼び出し前に評価されるため、片方に副作用や高コスト計算があると常に実行される（短絡評価ではない）。ドキュメント用コメントブロックはない。

---

## 可視化 (`visualization/`)

衛星軌道の 3D 表示・地上軌跡・日食（影）領域を描画し、アニメーション再生まで行う可視化ツール群。上位の `visualizeOrbit` / `visualizeTLE` がエントリポイントで、その内部で軌道伝播 (`propagateOrbit`)・座標変換 (`eci2latlon`)・太陽/影判定 (`calcSunShadow`)・地球描画 (`drawEarth3D`, `setupGroundTrack`) を組み合わせる。位置履歴 `rI` は行=時系列サンプル・列=成分 (Nx3, ECI, km) の横ベクトル規約で、`const`（`orbitConst()` が返す軌道定数：`GE` [km^3/s^2]、`RE` [km]、`RS` [km] 等）を全体で共有する。角度は基本的に rad、緯度経度プロット時のみ deg に変換する。

### 関数一覧
| 関数 | 概要 |
|---|---|
| `visualizeOrbit` | 軌道要素 or 位置速度履歴ファイルから単一衛星軌道を 3D + 地上軌跡で可視化 (スクリプト的エントリ) |
| `visualizeTLE` | TLE ファイル / `.mat` 履歴ファイル（複数可）から衛星軌道をまとめて可視化 (スクリプト的エントリ) |
| `propagateOrbit` | 軌道要素と時間配列から慣性系の位置・速度履歴を計算（ケプラー伝播） |
| `eci2latlon` | ECI 位置と GMST 列から測地緯度・経度を計算 |
| `calcSunShadow` | 太陽位置と日照/影係数（日食判定）を計算 |
| `drawEarth3D` | テクスチャ + 太陽方向影オーバーレイ付きの 3D 地球を描画し hgtransform ハンドルを返す |
| `drawEarth` | ワイヤフレーム球にテクスチャを貼った 3D 地球を GMST 回転付きで描画 |
| `setupGroundTrack` | 地上軌跡プロットの背景（地球テクスチャ・軸範囲・ラベル）を設定 |

### エントリポイント

#### `visualizeOrbit`
```matlab
visualizeOrbit(oe, options)
```
衛星の軌道運動を可視化する。左に 3D 軌道（慣性系）、右に地上軌跡を並べ、日照/影で色分けし、アニメーション再生する。`oe` を与えれば軌道要素からケプラー伝播、`rvFile` を与えれば `.mat` の位置速度（または軌道要素）履歴を読み込む。
- **入力**: `oe` (1x6, 混在単位) — 軌道要素 [a(km), e(-), i(rad), RAAN(rad), omega(rad), f(rad)]。省略時は ISS 風のデフォルト `[6798.137, 0.0001, 0.9006, 0.5236, 0, 0]`
- **オプション** (name-value): `nOrbits` (1x1, 既定 2) — 描画する周回数; `nPoints` (1x1, 既定 500) — サンプル点数; `animate` (1x1 logical, 既定 true); `showShadow` (1x1 logical, 既定 true) — 日照/影の色分けと影オーバーレイ; `showGroundTrack` (1x1 logical, 既定 true); `startDate` (1x6, 既定 `[2026 8 14 12 0 0]`) — エポック [年 月 日 時 分 秒]; `animSpeed` (1x1, 既定 50) — 再生速度（1/animSpeed 秒 pause）; `rvFile` (char, 既定 '') — 位置速度履歴 `.mat` パス
- **出力**: なし（figure を生成）
- **Note**: `rvFile` の `.mat` には `tSpan` が必須で、`rI` (Nx3) か `oeArray` (Nx6) のいずれかを含む。`jd0`（エポック JD）や `oeFlag`（oeArray 第6列種別: 1=真近点角=既定, 0=平均近点角）は省略可。`oe=[]` を渡し `rvFile` を指定するとファイル読み込みモードになる。ライブラリパスは関数内で自動 `addpath` される。
- *See also*: `visualizeTLE`, `propagateOrbit`, `calcSunShadow`, `eci2latlon`, `drawEarth3D`, `setupGroundTrack`

#### `visualizeTLE`
```matlab
visualizeTLE(varargin)
```
TLE ファイルまたは位置速度履歴 `.mat` ファイル（複数指定・混在可）から複数衛星の軌道を読み込み、同一図に色分けして可視化する。TLE は `readTLE` で軌道要素に変換し、`.mat` 履歴は共通時間軸へ spline 補間する。
- **入力** (可変長): 先頭からファイル名（拡張子 `.txt` / `.tle` / `.TLE` / `.mat`）を任意個。以降を name-value オプションとして解釈。ファイル無指定時は `sample_tle.txt` を使用
- **オプション** (name-value): `nOrbits` (既定 2); `nPoints` (既定 500); `animate` (既定 true); `showShadow` (既定 true); `animSpeed` (既定 50)
- **出力**: なし（figure を生成）
- **Note**: 複数衛星は共通の周回数・点数でサンプリングし、最長周期 `Tmax` を基準に時間軸を張る。エポック `jd0` は最初に見つかったファイル/TLE の値を採用し、無ければ `2026/8/14 12:00:00` を既定とする。地上軌跡は常に表示（`showGroundTrack` オプションは無い）。`.mat` の必須/省略変数は `visualizeOrbit` と同様（加えて `name` で衛星名を指定可）。オプション名は `nOrbits/nPoints/animate/showShadow/animSpeed` のいずれかが現れた時点でそれ以降を name-value と見なす。
- *See also*: `visualizeOrbit`, `readTLE`, `propagateOrbit`, `calcSunShadow`, `eci2latlon`, `drawEarth3D`, `setupGroundTrack`

### 計算ヘルパー

#### `propagateOrbit`
```matlab
[rI, vI] = propagateOrbit(oe, tSpan, mu)
```
軌道要素と時間配列からケプラー運動で慣性系の位置・速度履歴を計算する。初期真近点角を平均近点角に変換し、平均運動 n で各時刻の平均近点角を進めてから `oe2rv` で位置速度に変換する。
- **入力**: `oe` (1x6, 混在単位) — 軌道要素 [a(km), e, i(rad), RAAN(rad), omega(rad), f(rad)]（f=初期真近点角）; `tSpan` (1xN, s) — 時間配列; `mu` (1x1, km^3/s^2) — 地球重力定数（呼び出し側は `const.GE` を渡す）
- **出力**: `rI` (Nx3, km) — 慣性系位置履歴; `vI` (Nx3, km/s) — 慣性系速度履歴
- **Note**: 平均近点角は `mod(·, 2π)` で折り返すため、履歴を線に繋ぐと周回境界で見かけ上の不連続が入りうる（点描画では問題にならない）。

#### `eci2latlon`
```matlab
[lat, lon] = eci2latlon(rI, GMSTarray, const)
```
ECI 位置履歴を GMST で ECEF に回して測地緯度・経度に変換する。各サンプルごとに Z 軸 -GMST の DCM を適用し `ecef2LatLonH` で緯度経度を得る。
- **入力**: `rI` (Nx3, km) — 慣性系位置履歴; `GMSTarray` (Nx1, rad) — 各時刻の Greenwich Mean Sidereal Time; `const` — 軌道定数（WGS-84 楕円体パラメータに使用）
- **出力**: `lat` (Nx1, rad) — 測地緯度; `lon` (Nx1, rad) — 経度
- **Note**: `rI` と `GMSTarray` の行数（サンプル数）は一致している必要がある。

#### `calcSunShadow`
```matlab
[sunPos, nu] = calcSunShadow(jdArray, rI, const)
```
各時刻の太陽位置（VSOP）を計算し、地球による日食（影）判定を行う。`vsopConst` → `sun` で太陽位置、`shadow` で影係数を得る。
- **入力**: `jdArray` (Nx1 or 1xN, day) — ユリウス日配列; `rI` (Nx3, km) — 慣性系衛星位置; `const` — 軌道定数（`RS`, `RE` を使用）
- **出力**: `sunPos` (Nx3, km) — 慣性系の太陽位置; `nu` (Nx1) — 日照/影係数（1=日照、0=本影、半影・環状食では 0〜1 の中間値）
- **Note**: `nu` は厳密な 0/1 ではなく `shadow` が返す遮蔽係数（0〜1）。色分け等では `nu > 0.5` を日照とみなす運用。太陽位置と `rI` は同一単位（km）である必要がある。

### 描画ヘルパー

#### `drawEarth3D`
```matlab
hEarthTransform = drawEarth3D(GMST0, const, sunDir)
```
`naturalEarth.jpg` をテクスチャに貼った 3D 地球を描画し、GMST0 回転を適用した hgtransform ハンドルを返す（アニメーションではこのハンドルの Matrix を更新して自転させる）。`sunDir` を与えると慣性系に固定した半透明の夜側影オーバーレイと太陽方向矢印を追加する。
- **入力**: `GMST0` (1x1, rad) — 初期 GMST; `const` — 軌道定数（`RE` を使用）; `sunDir` (1x3, 既定 `[]`) — 太陽方向単位ベクトル。空なら影なし
- **出力**: `hEarthTransform` — 地球自転用の hgtransform ハンドル
- **Note**: テクスチャはライブラリ内の絶対パス `visualization/naturalEarth.jpg` から読み込む。影オーバーレイは地球半径の 1.002 倍の球で、太陽方向との内積で夜側の透過率を決める（慣性系固定なので地球本体だけが hgtransform で回る）。

#### `drawEarth`
```matlab
drawEarth(GMST, alpha, const)
```
ワイヤフレーム球に `earth.jpg` テクスチャを貼った 3D 地球を、GMST 分だけ Z 回転して描画する。背景を黒にし FaceAlpha を指定できる簡易版地球描画。
- **入力**: `GMST` (1x1, rad) — 地球自転角; `alpha` (1x1) — 地球表面の不透明度 (FaceAlpha, 0〜1); `const` — 軌道定数（`RE` を使用）
- **出力**: なし（現在の figure/axes に描画）
- **Note**: テクスチャはカレントディレクトリ相対の `earth.jpg` を読むため、`visualization/` を作業ディレクトリにするかパス上に置く必要がある（`drawEarth3D` は絶対パスで読むのと対照的）。
- *See also*: `orbitConst`

#### `setupGroundTrack`
```matlab
setupGroundTrack()
```
地上軌跡プロットの背景を設定する。`naturalEarth.jpg` を経度 -180〜180 deg・緯度 -90〜90 deg に貼り、軸範囲・アスペクト比・白い軸ラベル（経度/緯度 [deg]）を整える。
- **入力**: なし
- **出力**: なし（現在の axes を設定）
- **Note**: テクスチャはライブラリ内の絶対パス `visualization/naturalEarth.jpg` から読み込む。呼び出し後に緯度経度散布図を重ねる前提。

---

## サンプルスクリプト (`examples/`)

`examples/` はライブラリの使い方を示す実行可能なサンプル集で、軌道・姿勢のダイナミクス伝播、ライトカーブによる姿勢推定、SRP 計算、軌道決定、可視化などをテーマ別サブフォルダに分けて収録している。各サンプルは MATLAB でスクリプト (`main*.m` / `example*.m`) を開いてそのまま実行すればよい（ライブラリ本体の各ディレクトリに `addpath(genpath(...))` 済みであることが前提。`visualization/exampleOrbits.m` のように冒頭で自らパスを追加するものもある）。

### サブフォルダ別

- **`HCWeq/`** — Hill–Clohessy–Wiltshire 方程式による相対軌道運動の伝播例（`mainHCW.m`, `eomHCW.m`）。
- **`MEKF/`** — 剛体回転運動の真値生成と、Multiplicative Extended Kalman Filter による姿勢推定のデモ（`mainEuler.m`, `testMEKF.m`, `mekf.m`, `eulerEom.m`）。
- **`attitudeEstLightCurve/`** — 測光（ライトカーブ）データからの姿勢・軌道推定を各種非線形フィルタで実装した例。UKF・Lie 群上 UKF・Cubature KF・平方根 Cubature KF を比較（`mainLcUKFq.m`, `mainLcUKFqOrbit.m`, `mainLcUKFlie.m`, `mainLcCKFq.m`, `mainLcSRCKF.m`, `mainAttiOrbit.m` ほか）。
- **`attitudeOrbit/`** — 姿勢と軌道を連成させた運動の伝播例（EGM2008 重力・SRP を含む `mainAttiOrb.m`）と、セルフシャドウイング計算のテスト群（`selfShadowingTest/`）・軌道伝播テスト `JTOP/`。
- **`atttitudeControl/`** — 3-2-1 オイラー角で指定した目標姿勢への姿勢制御シミュレーション例（`mainAttiCtrl.m`, `eomAtti321.m`）。(注: フォルダ名は `atttitude` と t が重複した綴り)
- **`dualQuaternionMotion/`** — デュアルクォータニオンによる並進・回転連成運動学／力学の例と、SRP 込みのバリアント。分離解との比較も含む（`mainDQmotion.m`, `mainDQmotionSRP.m`, `dqEOM.m`, `decEOM.m`）。
- **`object/`** — `.obj` 衛星形状ファイルの読み込み・表示と、セルフシャドウイング関数の検証・MEX 化の例（`exampleReadingObj.m`, `exampleSelfShadowing.m`, `exampleMexSelfShadow.m`）。
- **`orbitDetermination/`** — 3 観測からの初期軌道決定手法の例。Gauss 法・Gibbs 法・double-r 法（`exampleGauss.m`, `exampleGibbs.m`, `exampleDoubleR.m`）。
- **`srp/`** — 太陽輻射圧の計算例。box-wing 形状に対する複数 SRP モデル（simple / AS / CT）の比較と、SRP のコンター可視化（`exampleSRP.m`, `exampleSRPcontour.m`）。
- **`visualization/`** — ISS・静止・Molniya・太陽同期など各種軌道タイプの可視化例（`exampleOrbits.m`）。

---

## 補助ディレクトリ（SPICE・演習・テスト・参考資料）

### SPICE (`SPICE/`)

NASA SPICE Toolkit（Mice）を MATLAB から利用するためのユーティリティ群。主なスクリプトは `loadSpiceKernel.m` で、`cspice_kclear` で既存カーネルをクリアしたうえで leapseconds（`naif0012.tls`）・地球 PCK・座標系フレーム定義などのカーネルを `cspice_furnsh` でまとめてロードする。

### 演習 Knock (`50attitudeKnock/`, `50orbitKnock/`)

学習用の「50本ノック」演習スクリプト集で、`50attitudeKnock/` が姿勢計算（`attiK_x_xx.m`：慣性テンソルの主慣性モーメント、回転表現など）、`50orbitKnock/` が軌道計算（`orbitK_x_xx.m`：オイラー角↔回転行列変換など）の解答コードを収める。各ファイルは問題文をコメントに埋め込んだ 1 問 1 スクリプト形式で、本体ライブラリの関数（例：`dcm1axis`）の使用例も兼ねる。

### テスト (`tests/`)

`matlab.unittest` ベースのユニットテストを置くディレクトリ。`runAllTests.m` がライブラリを path に追加してテストを一括実行し、結果をサマリ表で返す（非再帰実行なので `tests/manual/` は対象外）。テストファイル名から、SRP（`testSRP`）・座標系フレーム（`testFrames`）・軌道要素（`testOrbitElements`）・時刻変換（`testTime`）・姿勢（`testAttitude`）の各モジュールをカバーしていることが分かる。`tests/manual/` は自動判定に載せづらい計算を目視確認するためのプロット比較スクリプト群（`verifyMoonPosition`、`verifyEGM2008`、`verifyITRF2GCRF` など）。

### 参考資料 (`000refs/`)

ライブラリの背景となる参考資料（吉村の JGCD 2023 論文 PDF、補足メモ等）を置く保管用ディレクトリ。


---


## 付録: コードとドキュメントの不整合 (要確認)

初回の自動監査で、各関数のドキュメントコメント（`%[text]`）と実装（`function` シグネチャ・`arguments` ブロック）を突き合わせ、**148 件**の食い違いを検出しました（個人ファイルは除外済み）。**マニュアル本文の記述は実装（コード）側に合わせています。**

### 対応状況 (2026-07-07)

| 区分 | 件数 | 状態 |
|---|---|---|
| 実行時に必ず失敗する未定義呼び出し | 4 | ✅ コード修正済み |
| ドキュメントコメントのズレ（引数名・単位・出力・タイトル・壊れた See also・LaTeX ゴミ） | 約 94 ファイル | ✅ 実装に同期済み |
| コードロジックのバグ（数式・分岐・既定値） | 7 | ✅ 修正済み |
| 誤検出（doc は実装と一致 / `arguments` ブロックの記述） | 14 | — 対応不要 |

> 🧪 上記の全変更は **MATLAB R2026a で検証済み**: ユニットテスト `tests/runAllTests.m` が **66/66 pass**、加えて修正関数（`hms2deg`/`obliquity`、`normRow`、`calcAreaObj`〔正方形=1・長方形=2・台形=2.5・三角形=0.5〕、`day2s`、`axiQsol`、`gFigs`/`figs4Paper` の既定値）の数値スモークテストも pass。

### ✅ 修正済み: コードロジックバグ（2026-07-07）

コメントでは直らない実行コード側の不具合。実装を確認のうえ最小限で修正済み。

| ファイル | 内容 | 修正 |
|---|---|---|
| `conversion/hms2deg.m` | ~~時角(hour,min,sec)→度の換算式 out = hour.*15 + (min.*60+sec)./3600 の min/sec 係数が標準 15*(hour+min/60+sec/3600) と一致しない。意図が不明なため要確認（推測修正はせず）。~~ | **→ hour 項の誤った ×15 を除去し DMS→deg に統一（obliquity.m と整合。既存呼び出しは hour=0 のため数値不変）** |
| `lightcurves/visualizeBRDF.m:48` | ~~`[~,~,D(i,j)] = lcAS(sat,s,v)` が第3出力 cs を D に代入。lcAS の戻りは [sat,cd,cs,D] なので正しくは `[~,~,~,D(i,j)]`（直下のコメント lcCT 呼び出しは正しい形）。~~ | **→ lcAS の第4出力 D を取得するよう `[~,~,~,D(i,j)]` に修正** |
| `math/normRow.m` | ~~`if (Anorm>eps)` が複数行入力で列ベクトル比較になり、ゼロ行と非ゼロ行が混在すると全体が未正規化のまま返る。要素ごと判定＋行別分岐が必要。~~ | **→ 行ごとの正規化をベクトル化（ゼロ行は保持）。混在入力でも各行が正しく正規化される** |
| `object/selfShadow.m:44` | ~~三角/四角メッシュ判定を先頭面 `faces(1,4)` の NaN 有無だけで行っており、面ごとに三角・四角が混在するモデルで面 j を誤判定。`faces(j,4)` を見るべき。~~ | **→ 四角判定を面 j 自身の `faces(j,4)` に修正＋`size(faces,2)>=4` ガード追加** |
| `object/calcAreaObj.m:27,29` | ~~短絡 `\|\|` を意図する箇所で要素ごと `\|` を使用。加えて四角面の面積を単一外積ノルム（平行四辺形近似）で算出しており厳密でない。~~ | **→ `\|`→短絡 `\|\|`/`&&` に修正（Nx3 で index エラー回避）。四角面積を三角形2枚の和で厳密化（長方形では従来値と一致）** |
| `utility/gFigs.m` | ~~既定値 nDisp=1 を与える `arguments` ブロックがコメントアウトされているため、`gFigs()` 無引数呼び出しで `nDisp` 未定義→`if nDisp==0` でエラー。help は既定1と記載（doc は正、コード側の不備）。~~ | **→ 無引数呼び出し時に `nDisp=1` を補完** |
| `utility/figs4Paper.m` | ~~`if nargin<1, asis=0; elseif nargin<2, contentType='vector'; end` の既定値設定が誤り。0引数時に contentType 未定義のまま fig4Paper へ渡りエラー。`elseif` を独立した `if` にすべき。~~ | **→ 既定値設定を独立した `if` に分離（0引数で contentType 未定義になる不具合を解消）** |

### 初回監査ログ（参考）

<details><summary>検出された全 148 件を展開（重要度順。✅=修正済みの実行時バグ。doc のズレは 2026-07-07 に実装へ同期済み）</summary>


**🔴 High（25 件）**

| 状態 | ファイル | 内容 |
|---|---|---|
| ✅ 修正済 | `axiQsol.m` | ~~line 26 で存在しない関数 qMult_(4, 1, y, qIni) を呼んでいる。ライブラリには qMult (末尾アンダースコアなし、attitude/qMult.m) しかなく、qMult_ はどこにも定義されていない。そのまま実行すると未定義関数エラーで落ちる。正しくは qMult(4, 1, y, qIni)。~~<br>**→ qMult_ → qMult に修正 (2026-07-07)** |
| — | `calcRayIntersect.m` | ヘッダのドキュメントコメントが入力として satName および sat 構造体一式 (vertices/faces/area/pos/Ca/Cd/Cs/n/shadowFlag) を列挙しているが、実際の function シグネチャは calcRayIntersect(sun, nJ, vertJ, vertI) の 4 引数であり、ドキュメントと全く対応していない。arguments ブロック (sun(1,3), nJ(:,3), vertJ(:,3), vertI(:,3)) が正しい仕様。出力 flag の説明も ## outputs として書かれていない。 |
| — | `ctM2.m` | doc コメント (行4-5) は `v` と `sunB` を 1x3 vector と記すが、実装は nFacet x 3 行列を前提に動く: `h = h ./ vecnorm(h,2,2)` および `dot(sat.normal,h,2)` / `dot(v,h,2)` は行方向 (dim=2) 集約で、複数 facet の行列入力を想定している。呼び出し元 srpApproxCT2 も rRef を nFacet x 3 として渡す。doc の次元記載が実際の入力次元と不一致。 |
| — | `ctrlDq.m` | ドキュメントコメントは引数を 'w1, w2'（角速度）と 'v1, v2'（並進速度）でグループ化して列挙しているが、実際のシグネチャ function [dqa, dqb] = ctrlDq(dtp, dq1, dq2, w1, v1, w2, v2) では端点ごとに w1, v1, w2, v2 の交互順で並ぶ。呼び出し順を誤りやすい引数順序の記述不一致。 |
| — | `doubleR.m` | ドキュメントコメント(line 5)は入力に `rRange`（8次方程式の初期推定範囲）を挙げているが、実際のシグネチャ `doubleR(t, aziele, rObs, mu, RE)`(line 14)に `rRange` は存在しない。doc の入力リストと実装の引数個数・名前が食い違っている。 |
| — | `doubleR.m` | シグネチャの第5引数 `RE`（地球半径。line16 の収束判定 `TOL=1e-8*RE` と line27-28 の初期推定 `2.0*RE`,`2.01*RE` に使用）がヘッダコメントに一切記載されていない。未文書の必須入力。 |
| ✅ 修正済 | `doubleR.m` | ~~line24-25 で `day2s(t1-t2)` により日→秒へ換算しているが、`day2s` はライブラリ内に定義が存在しない（`conversion/` にあるのは逆変換の `s2day.m` のみ）。パス上に別途 `day2s` が無いと実行時に未定義関数エラーになる。~~<br>**→ conversion/day2s.m を新規作成（days→seconds） (2026-07-07)** |
| — | `egm2008.m` | The '## how to use' example block calls `[EGM.Cnm, EGM.Snm] = readEGM2008('EGM2008_to2190_TideFree.txt', EGM.GEODEG)` (filename first arg, two outputs), which no longer matches readEGM2008's current single-struct-in / single-struct-out signature. Following the example would pass a filename where the EGM struct is expected. |
| ✅ 修正済 | `geodeticIGRF.m` | ~~line 22 で igrfs(time, ...) を呼び出しているが、igrfs という関数はライブラリ内のどこにも存在しない (grep で定義なし)。この関数は現状そのままでは実行不能。~~<br>**→ 外部 igrfs 依存を撤去し内蔵 igrf12（igrfsyn12 mex）へ切替 (2026-07-07)** |
| — | `gprCov.m` | 出力の名前・意味がドキュメントと実装で不一致。タイトルは 'covariance matrix of GPR'、実装の出力は `sigPred`(予測分散, m×1)だが、doc の `## outputs` は `yPred`「predicted mean value ... mxn matrix」と予測平均を宣言している (gprMean からの copy-paste 誤り)。 |
| — | `gprCov.m` | 入力引数の個数不一致。doc の `## inputs` は `xAst, xTrain, yTrain, Kinv, hypPara` の5個を列挙するが、実装シグネチャ `gprCov(xAst, xTrain, Kinv, hypPara)` は4引数で `yTrain` を取らない (関数本体でも yTrain は未使用)。 |
| — | `gprMean.m` | ドキュメントコメントの入力欄と実装シグネチャが大きく不一致。doc は `xAst, xTrain, yTrain, Kinv, hypPara` の5引数を列挙するが、実装は `gprMean(xAst, xTrain, xMean, yTrain, yMean, L, hypPara)` の7引数。doc にある `Kinv` は実装に存在せず (Cholesky 因子 `L` に置き換わっている)、実装の `xMean` `yMean` `L` が doc から欠落している。 |
| — | `jacciaBowman.m` | doc コメント (line 4) は第2引数を `alp: right ascension`, geocentric longitude, rad と記載するが、実際のシグネチャ (line 19) は `lon` で、コード内 (line 69-70: XLON = lon; SAT(1) = mod(GWRAS + XLON, 2*pi)) では経度 (longitude) として使われる。引数名と意味が doc と不一致。 |
| — | `kernelGauss.m` | 出力の名前・意味・次元がすべて不一致。実装の出力は `k`(スカラーのカーネル値)だが、doc の `## outputs` は `yPred`「predicted mean value for test data, Mxn matrix」となっている (gprMean/gprCov からの copy-paste 誤り)。 |
| — | `lcSimple.m` | ドキュメントの ## outputs は4出力 (fObs, cd, cs, D) を列挙するが、実際のシグネチャ `[sat, cd, cs] = lcSimple(...)` は3出力で D が存在しない。D は他モデルからの流用と思われる陳腐化した記述。 |
| — | `math/associatedLegendre.m` | ヘッダーコメント `### 正規化されたAssociated Legendre多項式` は「正規化された」多項式を計算すると宣言しているが、実装 legendreRecursive は Condon-Shortley 位相 (-1)^m を含む非正規化の標準漸化式であり、正規化係数 (Schmidt 半正規化や完全正規化) は一切適用していない。タイトルが計算内容と意味的に食い違う。 |
| — | `oe2los.m` | 'See also roe2losApprox' references a function that does not exist anywhere in the library (find -iname roe2losApprox.m returns nothing). The actual approximate-LOS counterpart is named roe2mappedLOS. Broken See-also reference. |
| — | `qKine.m` | doc の入力欄が陳腐化。実シグネチャは (scalar, q, w) で w=角速度だが、doc は入力に 'R, rotation matrix, 3x3 matrix' を挙げ w(角速度)を全く記載していない。さらに出力を 'q: quaternions 1x4' としているが実際の出力は qKine(=クォータニオン時間微分 q̇)。q2DCM 系からのコピペ残骸。 |
| — | `rad2arcs.m` | タイトルが 'converting arcsecond to radian' となっており、実装 (arcs = angle .* 3600 .* 180 ./ pi、rad→arcsec) と逆。arcs2rad からの流用ミスで変換方向が誤っている。 |
| — | `readEGM2008.m` | doc/signature mismatch: header declares inputs `fName` (string), `deg`, `varagin` and outputs `Clm`, `Slm`, but signature is readEGM2008(EGM, deg, normalized) returning a single struct EGM with fields .Cnm/.Snm. The coefficient filename is hardcoded inside ('EGM2008_to2190_TideFree.txt'), not an argument. First arg is a struct to fill, third is a normalize flag. |
| — | `readTLE.m` | doc names the SPICE-selection argument `spiceFlag` with values on/off, but the actual third argument is `tool` (char) with allowed members {'yoshimuLibrary','SPICE','MATLAB'} (default 'yoshimuLibrary'). Argument name and value set both differ. |
| — | `srpCTinterp.m` | ドキュメントコメントは optional 入力として `NDF` と `nMC` を宣言しているが、実際の function シグネチャ `srpCTinterp(sat, sunB, d, const, correctionPara)` および arguments ブロックにこれらの引数は存在しない（受け取らない）。陳腐化した doc。 |
| ✅ 修正済 | `ukf.m` | ~~実装バグ: 31 行目 `xEst = sum(wim .* X, 1)` が function シグネチャ `ukf(f, h, Q, R, lam, P, x, y, wm, wc, wic)` に存在しない変数 `wim` を参照している。`wm` の誤記と思われ、そのままでは実行時に 'Undefined function or variable wim' でエラーになる。~~<br>**→ 未定義 wim → wm（平均重み）に修正 (2026-07-07)** |
| — | `zxz2q.m` | zyx2q.m と同じ問題。実シグネチャ (scalar, phi, theta, psi) に対し、doc に 'q: quaternions nx4' 'p: quaternions nx4' 'output: quaternions nx4' という qMult 由来の余分な入出力行が残っている。 |
| — | `zyx2q.m` | 入力の doc が誤り。実シグネチャは (scalar, phi, theta, psi)。phi/theta/psi は記載されているが、それに加えて 'q: quaternions nx4' 'p: quaternions nx4' 'output: quaternions nx4' という qMult 由来の入出力行がそのまま残っており、存在しない引数 q,p を入力であるかのように示している。 |

**🟡 Medium（48 件）**

| 状態 | ファイル | 内容 |
|---|---|---|
| — | `asg.m` | See also が存在しない関数を参照している: `sgMult` はライブラリ全体に存在しない（`sg` は存在）。 |
| — | `asg.m` | doc の定義式が異方性 SG ではなく等方 SG の式（μ·exp(λ(x^T p − 1))）をそのままコピーしている。実装は exp(−lam·VX² − mu·VY²)·max(0,VZ)·c で、x/y で帯域幅を分ける異方性形。式が実装と一致していない。 |
| — | `butcherTable.m` | line 14 の 'See also dqrk, dqgi.' の参照先 dqrk / dqgi がライブラリに存在しない（find で 0 件）。関連関数リンクが切れている。 |
| — | `calcLocalFrame.m` | ドキュメントは法線フィールドを sat.n と記載しているが、実装 (line 33, 40) が参照するのは sat.normal。引数名の相違。加えて Ca/Cd/Cs/area/pos/shadowFlag を入力として列挙するが実際に読むのは faces/vertices/normal のみ。 |
| — | `calcNormalObj.m` | ドキュメントは sat.n を入力フィールドのように列挙しているが、法線は本関数の出力 n であり、実際に読む入力は sat.vertices と sat.faces のみ。sat.n という入力は存在しない（出力/入力の混同）。 |
| — | `calcOrbitalState.m` | No documentation comment header at all (the only .m file in orbit/ lacking a %[text] title/description block); function calcOrbitalState(oe, mu) is undocumented despite being used by gve. |
| — | `dqConj.m` | 出力変数名が invDq（inverse を示唆）だが、実装 [dq(:,1) -dq(:,2:4) ...] はベクトル部の符号反転のみで共役（conjugate）を返す。逆元ではない。タイトル・日本語説明は正しく『共役』だが、出力名 invDq が誤解を招く（本物の逆元は dqInv）。 |
| — | `earthVSOP87.m` | doc describes the second input as `const: orbital constant`, but the signature's second argument is `earthVSOP` (the VSOP87 coefficient table produced by vsopConst), not the orbitConst struct. |
| — | `eulerEom.m` | ヘッダコメント (line 7) が 'used in eulerMain.m' と参照するが eulerMain.m はライブラリに存在しない。加えて Input/Output のドキュメント欄 (line 4-5) が空で、状態ベクトル x = [q(4); w(3)] や MOI の説明が全く記載されていない。 |
| — | `examples/attitudeEstLightCurve/mainLcUKFqOrbit.m, mainLcUKFlie.m` | タイトルコメントに 'not yet' と明記された未完成サンプル。ライブラリ利用者向けの動作サンプルとしては未完成である旨をマニュアルで明示するか、対象外とするのが望ましい。同様に orbitDetermination/exampleDoubleR.m も 'not yet'。 |
| — | `fig4Paper.m` | 4つの入力引数 asis/fig/nFig/contentType を持つが `## inputs` / `## outputs` セクションが一切なく、どの引数も説明されていない。さらに関数先頭のインラインコメント『第3引数はオプション（既定は "vector"）』は誤り。既定 "vector" を取るのは第4引数 contentType であり、第3引数は nFig。 |
| — | `figs4Paper.m` | nargin による既定値設定にバグ: `if nargin<1, asis=0; elseif nargin<2, contentType='vector'; end` の構造上、引数0個で呼ぶと contentType が未設定のまま fig4Paper に渡され実行時エラーになる（引数1個でのみ contentType が既定化される）。加えて `## inputs`/`## outputs` セクションが無い。 |
| — | `gFigs.m` | ヘルプ/arguments ブロックは `nDisp (1,1) ... = 1`（既定1）と宣言しているが、arguments ブロックは実装ではコメントアウトされているため既定値が効かず、引数なし `gFigs()` はヘルプの記載に反して未定義変数エラーになる。 |
| — | `gc2jd.m` | 入力引数の説明が実装と食い違う。実装シグネチャは6入力 (year, month, day, hour, minute, second) だが、doc コメント (行2-7) は year/month/hour/minute/second と出力 jd しか列挙せず、3番目の入力 `day` の説明が欠落している。 |
| — | `generateEulerAngleKinematics.m` | スクリプト(関数ではない、clc/clear を含む)。doc は '## output dxdt: euler angle kinematics' と宣言するが、実際に dxdt という変数は生成されず、出力は記号行列 B(および latex(B))。出力名が実装と不一致。 |
| — | `gprCov.m` | `See also gpCov.` が参照する `gpCov` はライブラリに存在しない (ライブラリ内の実関数は `gprCov`)。全ライブラリを検索しても `gpCov.m` は見つからず、正しくは `gprCov` のはず。 |
| — | `gprMean.m` | `See also gpCov.` が参照する `gpCov` はライブラリに存在しない (正しくは `gprCov` と思われる)。 |
| — | `hms2deg.m` | 実装 out = hour*15 + (min*60 + sec)/3600 は標準的な時角→度換算 (15*(hour + min/60 + sec/3600)) と分・秒項の係数が一致しない。doc は 'hour, min, and sec angles to deg' とだけ述べ算式の妥当性を示さないため、計算結果が想定と異なる可能性（数式バグの疑い、要確認）。 |
| — | `igrf12.m` | igrfsyn12 の mex ファイルへの依存 (line 18) だが、environment/ にはソース igrfsyn12.c/.h のみでコンパイル済み mex がない。mex 実体は examples/attitudeOrbit/igrfsyn12.mexmaci64 にあるため、そのパスを追加しないと呼び出せない。note に「use the mex file generated from igrfsyn12.c」とあるが生成物の所在は不明記。 |
| — | `jacciaBowman.m` | 核となる密度計算 JB2008 と補助関数 finddays / IERS / JPL_Eph_DE430 / timediff / iauCal2jd / iauGmst06 / Mjday_TDB がライブラリ内に存在せず (JB2008・finddays は grep で未検出)、外部ツールボックス (Mahooti HPOP 系) が path 上に必要。doc の references は `## NA` で依存が明示されていない。 |
| — | `jd2gc.m` | `See also jd2fyear.` (行15) の参照先 jd2fyear.m がライブラリ内のどこにも存在しない (find で全 yoshimuLibrary を走査済み)。壊れた関連関数参照。 |
| — | `kernelGauss.m` | `See also gpCov, kernelGaussMat.` のうち `gpCov` はライブラリに存在しない (正しくは `gprCov`)。`kernelGaussMat` は存在する。 |
| — | `kernelGaussMat.m` | `See also gpCov.` が参照する `gpCov` はライブラリに存在しない (正しくは `gprCov` と思われる)。 |
| — | `lcAS.m` | ドキュメントは出力を3つ (sat, cd, cs) しか記述しないが、シグネチャ `[sat, cd, cs, D] = lcAS(...)` は4出力で D (異方性分布項) が未記載。 |
| — | `lcCT.m` | ドキュメントの ## outputs 先頭が `fObs: total reflectance BRDF` だが、実際の第1出力は sat 構造体で、反射光量は sat.fObs フィールドに入る。lcSimple.m も同様に第1出力を fObs と誤記(実体は sat)。 |
| — | `math/associatedLegendre.m` | 冒頭コメント `P_l^m(x) を計算する複数の方法を実装` は複数手法を謳うが、本文では `方法1: 再帰関係式` のみが実装されており「方法2」以降は存在しない。陳腐化した/誤解を招くドキュメント。 |
| — | `math/normRow.m` | ゼロ割回避の分岐 `if (Anorm > eps)` は N>1 行のとき Anorm が列ベクトルになり、条件は全行が eps 超のときのみ真。ゼロ行と非ゼロ行が混在すると条件全体が偽となり行列全体が未正規化で返るため、doc の宣言『行列の各行を正規化』が多行かつゼロ行混在ケースで成立しない。 |
| — | `moonELP.m` | Docstring `## inputs` lists two separate coefficient arguments (`ELPcoeffA, ELPcoeffB: ELP coefficients`), but the actual signature is `[lon,lat,r] = moonELP(jd, ELP)` — a single struct `ELP` with fields `ELP.a` and `ELP.b`. Argument name and count mismatch between doc and implementation. |
| — | `oe2rv.m` | doc '## outputs'/'## inputs' and the commented arguments block name the gravitational-constant argument `GE`, but the function signature uses `mu` (oe2rv(oe, flag, mu)). Argument name mismatch (the revision note even says 'GEを直接引数として設定'). |
| — | `plotStd.m` | ドキュメントの `## inputs` は t / meanVal / stdVal の3入力しか記載していないが、実際の function シグネチャは `plotStd(t_, meanVal, stdVal, lineColor, regionColor)` の5入力。lineColor（既定 'r'）と regionColor（既定 'k'）が未記載。 |
| — | `qExp.m` | ドキュメントコメント (line 1-2) が 'この関数の簡単な概要です／この関数の詳細な説明です' というテンプレート雛形のまま放置され、実際のクォータニオン指数写像の説明になっていない。inputs/outputs/note/references セクションも未記載。 |
| — | `rad2arcs.m` | 出力のドキュメント記述が '`rad`: angle converted, arcs' だが、実際の function シグネチャの出力変数名は `arcs` (function arcs = rad2arcs(angle))。出力名が doc と不一致。 |
| — | `rodrigues2q.m` | 入出力名の取り違え。実シグネチャは (scalar, rod)→q だが、doc は入力欄に 'q: quaternions nx4' のみを挙げており、実際の入力 rod(Rodrigues パラメータ)が未記載。q は出力。 |
| — | `roe2mappedLOS.m` | Doc says azi is '-T軸方向からR軸への角度' (angle from -T axis toward R), but the code computes azi = atan2(rel(:,2), rel(:,1)), i.e. the angle measured from the R axis in the R-T plane. This matches oe2los.m's correct description ('R-T平面のR軸からの角度') and contradicts roe2mappedLOS's own doc. The azimuth-reference description is wrong. |
| — | `rv2oe.m` | doc header lists 8 separate quantities (a, e, inc, raan, w, nu, M, mu) as if outputs, but the function returns a single nx6 matrix oe = [a, e, inc, raan, w, nu]. M is computed then discarded, mu is an input not an output, and u/trueLon are commented out. |
| — | `sg.m` | See also が存在しない関数を 2 つ参照: `sgMult` も `sg3` もライブラリ全体に存在しない。 |
| — | `showSC.m` | ドキュメントは法線を sat.n と記載するが、実装 (line 32) は sat.normal を参照する。引数名の相違。また Normal='on' 時に描画に使う sat.uu がドキュメントに未記載、逆にドキュメントの sat.Ca/Cd/Cs は本関数で未使用。 |
| — | `srckfCorrGain.m` | 入力次元の記載ミス: doc は `X: 2n x n matrix` / `Y: 2n x m matrix` だが、SR-CKF は中心点付き 2n+1 点を用いる((2n+1)×n / (2n+1)×m が正)。同じ doc 内の `wc: 2n+1 x 1 vector` は正しいため矛盾している。 |
| — | `srckfCov.m` | 入力次元の記載ミス: doc は `X: 2n x n matrix` だが、srckfSigma が生成するのは (2n+1)×n。関数自体は nSigma=size(X,1) で追従するが doc が不整合。 |
| — | `srckfSigma.m` | 出力次元の記載ミス: doc は `Xout: 2n x n matrix` だが、実装(26-27 行)は `[x, x+xi, x-xi]` の中心点付き 2n+1 点を生成し、コメントにも '(2n+1) x n' とある。正しくは (2n+1)×n。 |
| — | `srpAS.m / srpASuni.m / srpCT.m / srpCTuni.m / srpCTinterp.m` | `See also srpLps, readSC, orbitConst.` の `srpLps` がライブラリ内に存在しない（.m ファイル・function 定義とも見つからず）。5 ファイルで参照先が壊れている。 |
| — | `srpApproxCT.m` | doc コメント (行4) は `thetaN ... 1x3 vector` と記すが、`arguments (Input)` ブロックは `thetaN (:,1)` (列ベクトル) を宣言し、実装も scalar 的に [0, sin(thetaN), cos(thetaN)] として使う。1x3 という次元記載は誤り (実体は角度スカラ/列ベクトル)。 |
| — | `srpApproxCT2.m` | doc コメント (行4) は `thetaN ... 1x3 vector` と記すが、`arguments (Input)` は `thetaN (:,1)` を宣言し、実装は facet ごとに thetaN(i) をスカラ参照する (nFacet 要素の列ベクトルが正しい)。1x3 という次元記載が誤り。 |
| — | `srpCT.m` | 既定 nMC の記載が不一致。arguments ブロックのコメントは `nMC = 10^3` だが、実際の実装 `if nargin < 6, nMC = 10^4` は 10^4 を設定する。ドキュメントの既定値が誤り。 |
| — | `srpSimple.m` | 出力 `srpCsOut` の doc が『total specular SRP force ... Nm, 1x3 matrix』と力(force)なのに単位を Nm(トルク)と記載。実装では rRef(鏡面反射方向の力成分)の総和であり単位は N。単位記載ミス。あわせて `srpCdOut` の説明が『force』なのに揃わず紛らわしい。 |
| — | `sun.m` | `See also sunLonLatR, earthVSOP, vsopConst.` references `earthVSOP`, but no `earthVSOP.m` exists in the library. The actual VSOP position function is `earthVSOP87` (orbit/earthVSOP87.m). Dead See-also reference. |
| — | `ukf.m` | ドキュメント(旧ヘッダコメント)が `R, observation noise: mx1 vector` と記載しているが、R は観測ノイズ共分散行列 (mxm) で、45 行目相当の Pyy+R で共分散として使われる。次元記載が誤り。 |
| — | `zyz2q.m` | 実シグネチャ (scalar, phi, theta, psi) に対し、doc の入力に 'q: quaternions nx4' 'p: quaternions nx4' という存在しない引数行が残っている(qMult 由来)。出力は別途 '## output q: nx4' があるので出力は妥当。 |

**⚪ Low（75 件）**

| 状態 | ファイル | 内容 |
|---|---|---|
| — | `50orbitKnock/` | Named '軌道計算50本ノック' but currently contains only a single solution script (orbitK_1_00.m), whereas 50attitudeKnock has 10; the orbit set is essentially unpopulated. |
| — | `README.md` | README の「ファイル構成」に exampleOrbits.m を挙げているが、ディレクトリに exampleOrbits.m は存在しない（.m ソースは calcSunShadow/drawEarth/drawEarth3D/eci2latlon/propagateOrbit/setupGroundTrack/visualizeOrbit/visualizeTLE の 8 個のみ）。README のリストが陳腐化している。 |
| — | `asg.m` | 引数 `mu` の説明の揺れ: doc 冒頭リストでは `mu` を 'bandwidth along y-axis'（帯域幅）としており実装と一致するが、等方 `sg` では `mu` が amplitude。役割が別なので混同注意。amplitude は別引数 `c`。 |
| — | `asg.m` | 出力名の不一致: doc は出力を `asg`（Nx1 vector）と記すが、function シグネチャの出力名は `aSG`。また実際の出力は入力配列と同形（NxM）で Nx1 固定ではない。 |
| — | `au2km.m` | See also の 'km2AU' は実ファイル km2au.m を指す。macOS 案件-insensitive FS / MATLAB の case-insensitive 呼び出しでは解決するが、大文字小文字が実ファイル名 (km2au) と不一致。 |
| — | `butcherTable.m` | line 6 の出力セクション見出しが '## otuputs' と綴り誤り（正しくは '## outputs'）。ドキュメント構造マーカーの誤りで、doc パーサが outputs セクションを認識できない可能性がある。 |
| — | `calcAreaObj.m` | 三角/四角判定に論理演算子 \| (elementwise OR) を短絡演算 \|\| の位置で使用しており (line 27, 29)、また四角面の面積を単一外積のノルム (平行四辺形近似) で算出している点はドキュメントに明記がない。意味に影響しうるが致命的でない実装上の注記。 |
| — | `ckfCorrGain.m` | doc は `xEst: state vector: nx1 vector` としているが、同じライブラリの ckfCov.m は `xEst: 1xn vector` と記載しており、xEst の向きの規約がファイル間で不統一(ckfCorrGain は xEst=xEst(:) で吸収するが ckfCov は 1xn 前提の xEst'*xEst を使う)。 |
| — | `ctM.m` | doc タイトル・数式に MATLAB Live Script の LaTeX エスケープ由来のゴミが残る: 行7 の `$M({\bf v})=\frac{G({\bf v})F({\bf v})}{4} $` は二重バックスラッシュでエスケープされており、そのままでは正しくレンダリングされない (ctM2.m 行7 も同様)。意味は G·F/4。 |
| — | `ctrlDq.m` | ドキュメントの ## inputs / ## outputs / ## input などの見出しが無く、引数説明が %[text] のフラットな行で書かれている（arguments ブロックも存在しない）。他関数と比べて型・次元の機械的宣言が欠けており、dq1/dq2/w/v の次元（1x3, 1x8 等）がコメントから読み取りづらい。 |
| — | `dcm1axisY.m` | 日本語1行説明が誤り: 'x軸周りの回転行列を計算' とあるが本関数は Y 軸まわり(dcm1axisX からのコピペ残骸)。タイトル行は正しく 'about Y axis'。 |
| — | `dcm1axisZ.m` | 日本語1行説明が誤り: 'x軸周りの回転行列を計算' とあるが本関数は Z 軸まわり(dcm1axisX からのコピペ残骸)。タイトル行は正しく 'about Z axis'。 |
| — | `dcm2q.m` | 本文コメントが陳腐化した実体無しファイルを参照: '連続性も確保したい場合は dcm2qC.mlx を使う' とあるが dcm2qC.mlx は存在せず、実体は dcm2qContinuous.m。 |
| — | `dcmI2RTN.m` | stale plain-comment header `% Input: raan, inc, w, f` while the actual 4th input parameter is named `nu` (dcmI2RTN(raan, inc, w, nu)). |
| — | `doubleR.m` | タイトル直後の入力説明で `rObs` を『3x3 matrix』としているが、変数名の意味（観測点位置ベクトル、行=時点）自体は正しい。ただし gauss.m では同じ量を `obsECI` と命名しており、姉妹関数間で観測点位置の引数名が `rObs`/`obsECI` と不統一。 |
| — | `dqMultMat.m` | ドキュメントは入力 dq の次元を明示せず、arguments ブロックは dq (:,8) と書くが、実装は内部で qMultMat と 4x4 zeros/固定 8x8 出力を構成しており単一の dual quaternion（1x8）を前提とする。バッチ（nx8）入力には非対応で、記載された可変行数の許容と実装の想定が食い違う。 |
| — | `drawShadowZones.m` | 他ファイルの規約（%[text] # タイトル / 日本語説明 / ## inputs / ## outputs / ## note 等）に従うドキュメントコメントブロックを一切持たず、function 先頭のインラインコメントのみ。ライブラリのドキュメント抽出（タイトル行・inputs/outputs）に乗らない。 |
| — | `eAnomaly.m` | `See also eAnomaly` is a self-reference (refers to the function itself instead of a related one, likely meant meanAnomaly). |
| — | `earthFullRotQ.m` | Title 'Earth rotation including precession and nutation' and signature are essentially identical to earthNutationPrecessionQ.m; earthFullRotQ's note says 'NA' whereas the twin function documents 'Polar motion is not included'. Duplicate/near-identical functions with inconsistent notes; unclear which is canonical. |
| — | `eulerEom.m` | line 2 のヘッダコメントに文字化け (mojibake) がある（'12�� 2, 2013' — 日本語の日付表記が壊れている）。 |
| — | `examples/attitudeEstLightCurve/mainLcUKFlie.m` | revision コメントに 'See also orbitConst, mainATOM.' とあるが mainATOM という名のスクリプトは examples 配下に存在しない (mainAttiOrbit / mainAttiOrb はある)。壊れた相互参照の疑い。 |
| — | `examples/attitudeEstLightCurve/mainLcUKFq.m ほか` | 複数の main スクリプト冒頭に 'cls' コマンドがある (clc/clear の後)。'cls' は標準 MATLAB コマンドではなくユーザ定義関数と思われ、パスに無いと実行時エラーになる。要確認。 |
| — | `examples/atttitudeControl` | サブフォルダ名の綴りミス: 'atttitudeControl' と t が3つ (正しくは 'attitudeControl')。ライブラリ本体の attitude 系ディレクトリと綴りが不一致で、addpath やドキュメント参照で混乱を招く可能性。 |
| — | `gFigs.m` | revisions ヘッダの日付は 20250623 だが、関数内ヘルプおよびインラインコメントの日付は 2025-07-25 で不一致。 |
| — | `gFigs.m` | ドキュメントコメントの `See also .` が空（参照先関数名が欠落したダングリング See also）。 |
| — | `geocentric2Geodetic.m` | `See also teme2Mod, mod2J2000, obliquity, nutation`: teme2Mod does not exist in the library, and the listed frame-transform functions are unrelated to this geocentric->geodetic coordinate conversion (likely copy-pasted See-also). |
| — | `gibbs.m` | ヘッダコメント line5 の日本語説明中に LaTeX エスケープ由来のゴミ `${\\rm km^3/s^2}$` が残っている（意味は「km^3/s^2」）。表示上のノイズ。 |
| — | `gmst.m` | `See also jd2GAST` references a function that does not exist anywhere in the library (likely renamed to gast). |
| — | `hms2deg.m` | 入力 hour / min / sec に単位・次元の記載がなく（doc は 'hour:' 'min:' 'sec.' のみ）、出力のみ 'out: angle, deg' と単位付き。入力側のドキュメントが不完全。 |
| — | `ifelse.m` | ドキュメントコメントブロック（タイトル・説明・inputs/outputs）を一切持たない。ライブラリのドキュメント規約に未準拠。 |
| — | `itrf2gcrf.m` | `See also gcrf2itrf, utc2tt`: gcrf2itrf does not exist in the library (only itrf2gcrf/qITRF2gcrf are present). Also the doc describes an input `leapJD` but the signature's second arg is the whole EOP struct. |
| — | `jacciaBowman.m` | doc の ## outputs (line 10-11) が `temp`: temperature の単位を空欄、`rho`: air density も単位未記載。他ラッパ (jr1971 は [K]/[kg/m^3] を明記) と比べ単位が欠落。 |
| — | `jd2gc.m` | doc コメント (行2-8) が入力 jd と出力の各成分を `## inputs` / `## outputs` の見出しなしで羅列しており、他ファイルの規約 (jd2mjd 等の ## inputs/## outputs) と不整合。意味は追えるが構造が欠けている。 |
| — | `jd2jdT.m` | タイトル直下の説明行が空 (行2) で、`## inputs` / `## outputs` セクションも無い。入力 jd・出力 T (ユリウス世紀数) の記載がドキュメントに一切ない。 |
| — | `kernelGaussMat.m` | 入力 `xTrain` の型記述が doc で「dxn vector」となっているが、実際は d 個のデータ点 × n 次元の行列 (vector ではなく matrix)。 |
| — | `lcAS.m` | See also に自分自身 lcAS を含む(lcSimple.m/lc.m/magInv.m の See also も同様に呼び出し元自身を列挙する冗長参照)。無害だが自己参照。 |
| — | `mag.m` | ドキュメントは `d` を `nFacet x 1` と記すが、d は観測者-対象間距離で全面積算後のフラックスに作用する量(面ごとの次元ではない)。加えて doc の第1引数説明 (`fObs:`) が空でシグネチャ引数名 fobs と表記が不一致。 |
| — | `math/README.md` | README は wrapPi を『角度を [-pi, pi] に正規化』(閉区間) と記すが、関数の doc/実装 (`mod(lambda+pi,2*pi)-pi`) は半開区間 [-pi, pi) で +pi は -pi に写る。区間表記が関数実体と不一致。 |
| — | `math/associatedLegendre.m` | legendreRecursive の入力チェック (m の範囲チェック `m<0\|\|m>l` の error、\|x\|>1 の warning) がすべてコメントアウトされており無効。コメントは入力条件を保証するように読めるが実際には検証されない。 |
| — | `math/wrapPi.m` | wrapPi は他ファイルの %[text] 形式ドキュメント規約と異なり、標準 MATLAB help 形式 (先頭 % コメント) で書かれている。整合性チェック上は正しく機能するが、ライブラリのドキュメント抽出規約 (%[text] # タイトル / ## inputs 等) から外れており自動マニュアル生成で拾い方が変わる。 |
| — | `mean2Osc.m` | doc '## inputs' lists both `f: true anomaly` and `M: mean anomaly`, but the signature has a single 6th positional argument M (mean2Osc(n, e, i, Ome, w, M, const)); the extra `f` entry does not correspond to any argument. |
| — | `meanAngle.m` | ドキュメントコメント(%[text] ブロック)が一切無い唯一の関数。タイトル・入出力・単位の説明が欠落しており、他ファイルの規約から外れている。 |
| — | `mod2J2000.m` | `See also teme2Mod, ...` references teme2Mod, which does not exist anywhere in the library. |
| — | `moonELP.m` | `See also sun.` is a weak/off-target reference for a Moon ELP position routine; more useful targets (moonLonLatR, readELP, moon) are the actual call chain. `sun` is unrelated to this function's role. Minor. |
| — | `oe2los.m` | The commented arguments block declares chiefQbi (1,4), but chiefQbi is passed to qRotation applied to the nx3 relative-position vector (rD-rC), which for time-series inputs (chiefOE/deputyOE are nx6) requires an nx4 quaternion. The (1,4) dimension annotation is inconsistent with the vectorized usage and too restrictive. |
| — | `oe2rv.m` | `See also roe2losApprox` references a function that does not exist anywhere in the library. |
| — | `precession.m` | `See also oe2roe, nutationQ, precession` includes precession itself (self-reference). Same self-referencing See-also is copy-pasted into precessionDCM.m and precessionQ.m. |
| — | `q2dcm.m` | 入力 doc 'q: 1x4 vector' は正しいが、arguments ブロックは q (:,4)(nx4 を許容)と宣言。実装本体は q(1)..q(4) のスカラー要素アクセスでバッチ非対応。arguments の (:,4) は誤解を招く(単一クォータニオン専用)。 |
| — | `qAxisAngle.m` | See also の参照先 'qCon' がライブラリに存在しない(qConj の誤り)。qCon.m はどこにも無い。 |
| — | `qConj.m` | doc 入力 'q: quaternions, 1x4 vector' だが arguments ブロック・実装は q (:,4)(nx4 バッチ対応)。次元記載が実装より狭い(1x4 と限定)。同様に qInv.m も 'q: 1x4 vector' と記載しつつ実装は nx4 対応。 |
| — | `qErr.m` | See also の参照先 'qCon' がライブラリに存在しない(qConj の誤り)。 |
| — | `qMult.m` | See also の参照先 'qCon' がライブラリに存在しない(qConj の誤り)。 |
| — | `qMultMat.m` | See also の参照先 'qCon' がライブラリに存在しない(qConj の誤り)。加えて doc の入力 'q: 1x4 matrix' は実装がスカラー要素アクセスでバッチ非対応なので妥当だが、'qMat: 4x4 matrix' が ## output 見出し無しで入力欄に混在している(軽微)。 |
| — | `qMultMat.m` | q2dcm.m と同種。arguments で q (:,4) と宣言するが、実装は q(1)..q(4) のスカラーアクセスで単一クォータニオン専用。バッチ入力は不可。 |
| — | `readSC.m` | ドキュメントの outputs 記述が sat の一部フィールド (vertices/normal/faces/area/pos/Ca/Cd/Cs/sunlitFlag) に限られ、実装が実際に生成する uu/uv/qlb/F0/kappa/nu/nv/mCT/fObs/MOI/m/force/torque が未記載。矛盾ではないが出力の網羅性欠如。 |
| — | `roe2mappedLOS.m` | Documentation comment lacks the '## inputs' / '## outputs' section headers used by the other files; the outputs azi/ele and input GE are interleaved among the input descriptions rather than grouped, making the input/output boundary unclear (cosmetic/structural, not a semantic error). |
| — | `sb.m` | ドキュメントコメントの `See also .` が空（参照先関数名が欠落したダングリング See also）。 |
| — | `sclerp.m` | ドキュメント本文に ## note / ## references / ## revisions / See also の各セクションが欠落しており、他ファイルのテンプレート構造と不整合（内容の誤りではなくドキュメント構造の欠落）。 |
| — | `selfShadow.m` | See also selfShadow がヘッダの自己参照になっている（自分自身を See also に列挙）。文脈上、下請けの calcRayIntersect を指すのが妥当。 |
| — | `selfShadow.m` | 四角/三角メッシュの分岐を先頭面の 4 列目 faces(1,4) の NaN 有無だけで判定しており (line 44)、面ごとに三角・四角が混在するモデルでは面 j が四角でも三角として扱われる/その逆が起こりうる。calcAreaObj が行ごとに判定しているのと不整合。 |
| — | `sg.m` | 出力次元の記載ミス: doc は `G` を 'Nx1 vector' とするが、実装 `mu.*exp(lam.*(...))` は入力 x,y,z（NxM 行列）と同形の配列を返す。 |
| — | `sgMix.m` | 次元記載の不一致: doc の inputs&outputs は `lam*`/`mu*` を 'scalar' と書くが、arguments ブロックは `lam1 (:,1)`/`mu1 (:,1)` の縦ベクトル（n 個一括処理可）を宣言し、実装も vecnorm(...,2,2) で nx1 を返す。実際にはスカラーに限定されない。 |
| — | `srpApproxCT.m` | `See also srpApproxCT.` (行15) が自分自身を参照している (自己参照)。srpApproxCT2 または ctM/ctM2 を指すべき See also の写し間違いと思われる。 |
| — | `srpApproxCT.m` | doc コメント (行5) は `sunB ... 1x3 vector` と記すが、`arguments (Input)` は `sunB (:,3)` を宣言する (複数行許容)。単一 facet 用途では 1x3 で整合するが、次元記載の厳密性としては arguments と不一致。 |
| — | `srpApproxCT2.m` | `See also srpApproxCT.` (行15) の下、doc タイトル・本文は srpApproxCT からのコピーで、この関数が複数 facet 版である旨や sunlitFlag が 1 固定 (行42、日陰判定無効) である旨がドキュメントに反映されていない。srpApproxCT (行40) は太陽との内積で sunlitFlag を計算しており、2 関数の挙動差が doc に記載されていない。 |
| — | `srpCTinterp.m` | 出力 `srpCdOut`/`srpCsOut` を doc は『total diffuse/specular part of SRP』と単純総和のように記すが、実装は総力から拡散・入射(impinged)成分を差し引いた近似分離（コメントにも Approximate と明記）。他の srp* 関数の厳密な成分分離とは意味が異なる。 |
| — | `srpCTuni.m` | input ドキュメントで第5引数を ``varagin``（可変長引数）と呼び『NDF distribution function, default: Beckmann』と説明するが、実際のシグネチャは名前付き固定引数 `NDF`（+ `nMC`）。引数名・種別の記載が実装と不一致。 |
| — | `sunG.m` | Docstring input header uses `## input` (and output uses ``## `output` `` with backticks) whereas the codebase convention/prompt spec is `## inputs` / `## outputs`. Cosmetic doc-format drift, not a semantic error. |
| — | `sunLonLatR.m` | Output `r` is in AU whereas the analogous distance output of the Moon routines (moonLonLatR/moonELP `r`) and all position vectors are in km. The doc does correctly state AU, but the cross-function unit inconsistency is an easy misuse trap; `sun.m` relies on `au2km` to compensate. |
| — | `teme2J2000.m` | `See also teme2Mod, ...` references teme2Mod, which does not exist in the library. Also the doc '## inputs' section omits `const`, but the signature requires const as the 5th argument. |
| — | `tod2Mod.m` | `See also teme2Mod, mod2J2000, obliquity, nutation`: teme2Mod does not exist in the library. |
| — | `ukf.m` | ヘッダコメント 'related function files: ukf_sigma, ukf_cov, ukf_corr, testUKF' が参照する 4 関数はいずれもライブラリに存在しない(現行の ukfSigma / ukfCov / ukfCorrGain へ改名済み)。陳腐化した See also 参照。 |
| — | `ukfCorrGain.m` | ドキュメントコメントに `## input` / `## output` の区切り見出しが欠落しており、入力(xEst,X,yEst,Y,wc,R)と出力(Pyy,Pxy)が同一ブロックに混在。さらに 3 つ目の出力 K がコメントに記載されていない(シグネチャは [Pyy, Pxy, K])。 |
| — | `utc2tt.m` | `See also orbitConst` (行14) の参照先は orbit/orbitConst.m に実在するが、これは軌道定数を返す関数で時刻変換とは無関係。jd2mjd.m / mjd2jd.m / ut2tt.m も同じ orbitConst を See also に挙げており、関連度の低い定型参照が使い回されている疑い。 |
| — | `visualizeBRDF.m` | line 48 `[~, ~, D(i,j)] = lcAS(sat, s, v)` は lcAS の第3出力 cs (鏡面成分) を D に代入している。lcAS の分布項 D は第4出力なので、可視化される量は意図した分布項ではなく鏡面反射成分。直下のコメントアウト lcCT 呼び出し (line 49) は正しく `[~,~,~,D(i,j)]` としており不整合。 |

</details>

