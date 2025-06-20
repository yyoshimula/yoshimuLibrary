%[text] ## 軌道計算50本ノック
%[text] ## 1-00
%[text] 適当な3-1-3オイラー角を設定し，回転行列へ変換せよ．
clc
clear
close all
%%
%[text] ## 各軸周りの回転行列
%[text] ここでは無名関数使ってみる．普通に関数を作るのでもok．
%[text] cf. yoshimuLibraryのdcm1axis.mlx
% x-axis
Rx = @(phi)[1 0 0
    0 cos(phi) sin(phi)
    0 -sin(phi) cos(phi)];

% y-axis（使わないけど．．）
Ry = @(phi)[cos(phi) 0 -sin(phi)
    0 1 0
    sin(phi) 0 cos(phi)];

% z-axis
Rz = @(phi)[cos(phi) sin(phi) 0
    -sin(phi) cos(phi) 0
    0 0 1];

%%
%[text] ## 回転行列
%[text] 回転角は$\\phi \\rightarrow \\theta \\rightarrow\\psi$の順番とする．
phi = 0.1 %[output:661d9a1c]
theta = 0.3;
psi = -0.3;

R = Rz(psi) * Rx(theta) * Rz(phi) %[output:77d2cac4]


%[text] ## おまけ
%[text] cf zxz2dcm.mlx
%[text] ↑と比較して検算
R2 = zxz2dcm(phi, theta, psi) %[output:51a2e2a5]
%[text] 回転行列は$R^{-1} = R^T\n$なので，
R * R2' %[output:19a412a7]
%[text] を計算して単位行列になるかを確認するのでも良い．

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:661d9a1c]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"phi","rows":"1","value":"0.1000"},"version":0}
%---
%[output:77d2cac4]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"R","rows":3,"type":"double","value":[["0.9787","-0.1855","-0.0873"],["0.2029","0.9376","0.2823"],["0.0295","-0.2940","0.9553"]]}}
%---
%[output:51a2e2a5]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"R2","rows":3,"type":"double","value":[["0.9787","-0.1855","-0.0873"],["0.2029","0.9376","0.2823"],["0.0295","-0.2940","0.9553"]]}}
%---
%[output:19a412a7]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":3,"type":"double","value":[["1.0000","-0.0000","-0.0000"],["-0.0000","1.0000","0.0000"],["-0.0000","0.0000","1.0000"]]}}
%---
