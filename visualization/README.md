# OLC - Orbit Visualization for MATLAB

衛星の軌道運動を可視化するMATLABプログラム集

## 概要

yoshimuLibraryの軌道計算・可視化関数を活用した衛星軌道シミュレータです。
3D軌道表示、地上軌跡、日食領域の可視化、アニメーション機能を備えています。

## 必要環境

- MATLAB R2019b以降（推奨）
- yoshimuLibrary（`/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary`）
- Mapping Toolbox（地上軌跡表示用、オプション）

## ファイル構成

```
OLC/
├── README.md           # このファイル
├── visualizeOrbit.m    # メイン可視化関数
├── exampleOrbits.m     # 様々な軌道タイプの例
├── visualizeTLE.m      # TLEファイルからの可視化
└── sample_tle.txt      # サンプルTLEデータ
```

## クイックスタート

```matlab
% ライブラリパスを追加（自動で行われます）
addpath(genpath('/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary'));

% ISS風の軌道をデフォルトで表示
visualizeOrbit()
```

## 使い方

### 1. 基本的な使い方（visualizeOrbit）

```matlab
% デフォルト軌道（ISS: 高度420km, 傾斜角51.6度）
visualizeOrbit()

% カスタム軌道要素を指定
% oe = [a, e, i, RAAN, omega, f]
%   a: 軌道長半径 [km]
%   e: 離心率 [-]
%   i: 軌道傾斜角 [rad]
%   RAAN: 昇交点赤経 [rad]
%   omega: 近点引数 [rad]
%   f: 真近点角 [rad]

const = orbitConst();
oe = [const.RE + 800, 0.001, deg2rad(98), deg2rad(45), 0, 0];
visualizeOrbit(oe)
```

### 2. ファイルからの読み込み（visualizeOrbit）

位置速度履歴または軌道要素履歴を`.mat`ファイルで渡して可視化できます。

```matlab
% 位置速度履歴ファイルから読み込み
visualizeOrbit([], rvFile='data.mat')

% 開始日時を指定して読み込み
visualizeOrbit([], rvFile='data.mat', startDate=[2026 3 14 12 0 0])
```

`.mat`ファイルには以下の変数を含めてください（`rI`または`oeArray`のいずれか必須）:

| 変数 | 必須 | 説明 |
|------|------|------|
| `rI` | △ | 位置履歴 [km], Nx3 (ECI) |
| `oeArray` | △ | 軌道要素履歴 [a,e,i,RAAN,omega,f], Nx6 |
| `tSpan` | 必須 | 時刻配列 [s], 1xN or Nx1 |
| `jd0` | 任意 | エポックのユリウス日（省略時は`startDate`を使用） |
| `oeFlag` | 任意 | `oeArray`第6列の種別（1=真近点角(default), 0=平均近点角） |
| `vI` | 任意 | 速度履歴 [km/s], Nx3 |

### 3. オプション指定

```matlab
visualizeOrbit(oe, ...
    'nOrbits', 3, ...           % 表示する周回数
    'nPoints', 1000, ...        % プロット点数
    'animate', true, ...        % アニメーションON/OFF
    'showShadow', true, ...     % 日食領域の表示
    'showGroundTrack', true, ...% 地上軌跡の表示
    'animSpeed', 100, ...       % アニメーション速度
    'startDate', [2026 3 14 12 0 0])  % 開始日時
```

### 4. 様々な軌道タイプ（exampleOrbits）

```matlab
exampleOrbits
```

対話形式で以下の軌道を選択できます:
1. **ISS** - 低軌道（高度420km, 傾斜角51.6度）
2. **GEO** - 静止軌道（高度35,786km）
3. **Molniya** - 高離心率楕円軌道（周期12時間, 傾斜角63.4度）
4. **SSO** - 太陽同期軌道（高度700km, 傾斜角98度）
5. **GPS** - 中軌道（周期約12時間, 傾斜角55度）
6. **カスタム** - 自由に軌道要素を入力
7. **複数軌道同時表示** - 上記軌道を一度に表示

### 5. TLEからの可視化（visualizeTLE）

```matlab
% サンプルTLEファイルを使用
visualizeTLE('sample_tle.txt')

% 複数ファイル
visualizeTLE('iss.tle', 'hubble.tle')

% オプション指定
visualizeTLE('sample_tle.txt', 'nOrbits', 3, 'animate', true)
```

### 6. ファイルからの読み込み（visualizeTLE）

`.mat`ファイルをTLEファイルと同様に引数で渡せます。TLEと`.mat`の混在も可能です。

```matlab
% 位置速度履歴ファイルから可視化
visualizeTLE('orbit.mat')

% TLEファイルと履歴ファイルを同時に表示
visualizeTLE('iss.tle', 'orbit.mat')

% 軌道要素履歴ファイルも対応
visualizeTLE('oe_history.mat', 'nOrbits', 3)
```

`.mat`ファイルの仕様は`visualizeOrbit`と共通です。追加で以下の変数が使えます:

| 変数 | 必須 | 説明 |
|------|------|------|
| `name` | 任意 | 衛星名（省略時はファイル名を使用） |

## 機能詳細

### 3D軌道表示
- 慣性座標系（J2000）での軌道表示
- 地球テクスチャマッピング
- 太陽方向の矢印表示
- 日照/影領域の色分け（黄色: 日照、青紫: 影）

### 地上軌跡（Ground Track）
- ECEF座標系での軌跡表示
- 海岸線の重ね合わせ
- 緯度経度グリッド

### アニメーション
- リアルタイム衛星位置追跡
- 時刻表示
- 速度調整可能

### 日食計算
- 本影（umbra）・半影（penumbra）の判定
- VSOP87による高精度太陽位置計算

## 使用ライブラリ関数

| カテゴリ | 関数 | 説明 |
|---------|------|------|
| 軌道計算 | `oe2rv()` | 軌道要素→位置・速度変換 |
| | `orbitConst()` | 軌道定数（地球半径、重力定数等） |
| | `readTLE()` | TLEファイル読み込み |
| 可視化 | `drawEarth()` | 地球テクスチャ表示 |
| 座標変換 | `gmst()` | グリニッジ恒星時 |
| | `qRotation()` | クォータニオンによる回転 |
| 太陽・月 | `sun()` | 太陽位置（J2000座標系） |
| | `vsopConst()` | VSOP87係数 |
| 影計算 | `shadow()` | 日食判定 |
| 時刻 | `gc2jd()` | グレゴリオ暦→ユリウス日 |
| | `jd2gc()` | ユリウス日→グレゴリオ暦 |

## 軌道要素の説明

```
        ↑ Z (北極方向)
        |
        |    軌道面
        |   /
        |  /  ω (近点引数)
        | /______ 近地点
        |/     ↗
   -----+--------→ X (春分点方向)
       /|\  Ω (昇交点赤経)
      / | \
     /  |  \
    /   |   \ i (軌道傾斜角)
   Y    |
        ↓ 昇交点
```

| 要素 | 記号 | 単位 | 説明 |
|------|------|------|------|
| 軌道長半径 | a | km | 楕円の長径の半分 |
| 離心率 | e | - | 楕円の扁平度（0=円, 1=放物線） |
| 軌道傾斜角 | i | rad | 赤道面との角度 |
| 昇交点赤経 | Ω | rad | 春分点から昇交点への角度 |
| 近点引数 | ω | rad | 昇交点から近地点への角度 |
| 真近点角 | f | rad | 近地点からの衛星位置角度 |

## トラブルシューティング

### 地球テクスチャが表示されない
`earth.jpg`がMATLABパス上にあることを確認してください。
yoshimuLibrary/object/フォルダ内にあります。

### coastlinesが見つからない
Mapping Toolboxがインストールされていない場合、地上軌跡の海岸線は表示されません。
軌跡自体は問題なく表示されます。

### アニメーションが遅い
`animSpeed`パラメータを大きくするか、`nPoints`を減らしてください。

## ライセンス

yoshimuLibraryと同様のライセンスに従います。

## 作成日

2026-03-14
