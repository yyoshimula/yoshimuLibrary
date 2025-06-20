%[text] ## 姿勢計算50本ノック
%[text] ## 2-00
%[text] 1-00で設定した慣性テンソルJの主慣性軸への回転行列に対応する1-2-3オイラー角を求めよ
clc
clear
close all
%%
%[text] ## 慣性テンソル
%[text] 対称行列であることに注意
J = [10 2 -3
    2 14 1
    -3 1 12];
%%
%[text] ## 主慣性モーメント
[V, D] = eig(J); % 固有ベクトルV, 固有値Dを求めるeig関数．

%[text] 固有ベクトルを並べたものは慣性テンソルを次式で対角化する．
%[text] テンソルの座標変換でもあるので，$V\n$そのものが回転行列
V^(-1) * J * V %[output:26ca54fb]

%[text] ## 1-2-3オイラー角
%[text] $\\phi \\rightarrow \\theta \\rightarrow \\psi$ の回転順とする．回転行列は
%[text] $R = \\left\[\\begin{array}{ccc}\\cos(\\theta)\\cos(\\psi) & \\sin(\\theta)\\cos(\\psi)\\sin(\\phi)+\\sin(\\psi)\\cos(\\phi) &\\sin(\\psi)\\sin(\\phi)-\\sin(\\theta)\\cos(\\psi)\\cos(\\phi) \\\\\n-\\cos(\\theta)\\sin(\\psi) &\\cos(\\psi)\\cos(\\phi)-\\sin(\\theta)\\sin(\\psi)\\sin(\\phi) &  \\sin(\\theta)\\sin(\\psi)\\cos(\\phi)+\\cos(\\psi)\\sin(\\phi)\\\\\n  \\sin(\\theta) &-\\cos(\\theta)\\sin(\\phi) & \\cos(\\theta)\\cos(\\phi) \n\\end{array}\n\\right\]$
%[text] なので，$R\_{3,1}, R\_{3,2}, R\_{3,3}, R\_{1,1},R\_{2,1}$から逆算

R = V;

phi = atan2(-R(3,2), R(3,3)) % rad %[output:2f86622e]
theta = atan2(R(3,1), sqrt(R(1,1).^2 + R(2,1).^2)) %[output:5f73cbc2]
psi = atan2(-R(2,1), R(1,1)) %[output:181f5391]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:26ca54fb]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":3,"type":"double","value":[["7.1270","-0.0000","-0.0000"],["-0.0000","14.0000","0"],["0.0000","-0.0000","14.8730"]]}}
%---
%[output:2f86622e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"phi","rows":"1","value":"-1.8039"},"version":0}
%---
%[output:5f73cbc2]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"theta","rows":"1","value":"0.5750"},"version":0}
%---
%[output:181f5391]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"psi","rows":"1","value":"0.3738"},"version":0}
%---
