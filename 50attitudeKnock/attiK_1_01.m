%[text] ## 姿勢計算50本ノック
%[text] ## 1-01
%[text] 適当な2次元回転行列を設定し，適当な2次元ベクトルaの回転を計算せよ．
clc
clear
close all

%[text] ## 2次元ベクトル $\\bf a\n$
a = [0
    2]; % 横ベクトル（row vector）でもいい
a = a ./ norm(a); % 正規化（しなくてもいい）

%[text] ## 回転行列
%[text] 問題文はあえて曖昧に書いている．
%[text] 座標系の回転とベクトルの回転の区別を理解した上で計算しているならどちらでもok（この意味がよくわかんない場合は[処方箋](https://www.overleaf.com/read/fsxbnvdjjhdp)を読む）
%[text] $R\_c = R^T\_v$の関係性がある．
phi = 0.3; % rad, 回転角

%[text] ### 回転行列（座標系の回転）
%[text] $R\_c = \\left\[ \\begin{array}{cc}\n\\cos{\\phi} & \\sin{\\phi} \\\\\n-\\sin{\\phi} & \\cos{\\phi}\n\\end{array} \\right\]$
Rc = [cos(phi) sin(phi)
    -sin(phi) cos(phi)];
b = Rc * a 
%[text] ### 回転行列（ベクトルの回転）
%[text] $R\_v = \\left\[ \\begin{array}{cc}\n\\cos{\\phi} & -\\sin{\\phi} \\\\\n\\sin{\\phi} & \\cos{\\phi}\n\\end{array} \\right\]$
Rv = [cos(phi) -sin(phi)
    sin(phi) cos(phi)];

b = Rv * a
b = Rc' * a

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
