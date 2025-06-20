%[text] ## 姿勢計算50本ノック
%[text] ## 2-01
%[text] 2-00で求めた1-2-3オイラー角を3-2-1オイラー角へ変換せよ
clc
clear
close all
%%
%[text] ## 1-2-3オイラー角
%[text] $\\phi \\rightarrow \\theta \\rightarrow\\psi\n$ の順とする．ここでは適当に値を設定
phi = 0.1;
theta = 0.2;
psi = 0.3;
R = [cos(theta)*cos(psi) sin(theta)*cos(psi)*sin(phi)+sin(psi)*cos(phi) sin(psi)*sin(phi)-sin(theta)*cos(psi)*cos(phi)
    -cos(theta)*sin(psi) cos(psi)*cos(phi)-sin(theta)*sin(psi)*sin(phi) sin(theta)*sin(psi)*cos(phi)+cos(psi)*sin(phi)
    sin(theta) -cos(theta)*sin(phi) cos(theta)*cos(phi)];
%%
%[text] ## 3-2-1オイラー角の回転行列
%[text] $\\phi \\rightarrow \\theta \\rightarrow\\psi\n$ の順でz軸まわり → y軸まわり → x軸まわりの順とする．
%[text]  資料によっては$\\phi\n$がx軸周りだったり，定義が異なる場合があるので注意．
%[text] $R=\\left\[\\begin{array}{ccc}\n\\cos{\\theta}\\cos{\\phi} & \\cos{\\theta}\\sin{\\phi} & -\\sin{\\theta}\\\\\n\\sin{\\theta}\\sin{\\psi}\\cos{\\phi}-\\cos{\\psi}\\sin{\\phi} & \\sin{\\theta}\\sin{\\psi}\\sin{\\phi}+\\cos{\\psi}\\cos{\\phi} & \\cos{\\theta}\\sin{\\psi}\\\\\n\\sin{\\theta}\\cos{\\psi}\\cos{\\phi}+\\sin{\\psi}\\sin{\\phi} & \\sin{\\theta}\\cos{\\psi}\\sin{\\phi}-\\sin{\\psi}\\cos{\\phi} & \\cos{\\theta}\\cos{\\psi}\n\\end{array}\\right\]$
%[text] 回転行列はどのような姿勢表現パラメータを使っても一意なので，1-2-3オイラー角から回転行列を計算し，$R\_{3,2},R\_{3,3},R\_{1,1},R\_{1,2}$から3-2-1オイラー角を逆算
psi = atan2(R(2,3), R(3,3)) %[output:5614d3ee]
phi = atan2(R(1,2), R(1,1)) %[output:1940a75e]
theta = atan2(-R(1,3), sqrt(R(2,3).^2 + R(3,3).^2)) %[output:6e540f9b]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:5614d3ee]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"psi","rows":"1","value":"0.1564"},"version":0}
%---
%[output:1940a75e]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"phi","rows":"1","value":"0.3226"},"version":0}
%---
%[output:6e540f9b]
%   data: {"dataType":"not_yet_implemented_variable","outputData":{"columns":"1","name":"theta","rows":"1","value":"0.1600"},"version":0}
%---
