%[text] ## 姿勢計算50本ノック
%[text] ## 1-02
%[text] 適当な3次元回転行列を設定し，適当な3次元ベクトルbの回転を計算せよ．
clc
clear
close all

%[text] ## 3次元ベクトル
b = [1
    2
    3]; % 横ベクトル（row vector）でもいい
b = b ./ norm(b); % 正規化（しなくてもいい）
%%
%[text] ## 回転行列 R
%[text] 問題文はあえて曖昧に書いている． 座標系の回転とベクトルの回転の区別がわかって計算しているならどちらでもok
phi = 0.1; % rad, 回転角
theta = 0.2; % rad, 回転角
psi = 0.3; % rad, 回転角

%[text] ### 座標系の回転
% x軸周り
Rx = [1 0 0
    0 cos(phi) sin(phi)
    0 -sin(phi) cos(phi)];

% y軸周り
Ry = [cos(theta) 0 -sin(theta)
    0 1 0
    sin(theta) 0 cos(theta)];

% z軸周り
Rz = [cos(psi) sin(psi) 0
    -sin(psi) cos(psi) 0 
    0 0 1];

%[text] ３次元回転行列であればどういうものでもok
%[text] 回転行列を1つだけ使うのでもok
c = Rx * Ry * Rz * b
%[text] ### ベクトルの回転
%[text] 省略
%%
%[text] ## おまけ 無名関数を使ってみる
%[text] たまに簡単な関数を定義したいときは無名関数を使う→ いちいち関数ファイルを作るまでもない時に便利
DCMx = @(x) [1 0 0
    0 cos(x) sin(x)
    0 -sin(x) cos(x)];

DCMy = @(y)  [cos(y) 0 -sin(y)
    0 1 0
    sin(y) 0 cos(y)];

DCMz = @(z) [cos(z) sin(z) 0
    -sin(z) cos(z) 0 
    0 0 1];

c = DCMx(phi) * DCMy(theta) * DCMz(psi) * b

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
