%[text] ## 姿勢計算50本ノック
%[text] ## 2-02
%[text] 2-01で求めた3-2-1オイラー角をクォータニオンへ変換せよ
clc
clear
close all
%%
%[text] ## 3-2-1オイラー角の回転行列
%[text] $\\phi \\rightarrow \\theta \\rightarrow\\psi\n$ の順とする．ここでは適当に値を設定
%[text] 資料によっては$\\phi$がx軸周りだったり，定義が異なる場合があるので注意．
phi = 0.1;
theta = 0.2;
psi = 0.3;

R = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
     sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(psi)
     sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) cos(theta)*cos(psi)];
%%
%[text] ## 回転行列から逆算する方法（もあるが面倒）
%[text] 同様にクォータニオンの回転行列が
%[text] $R\_{b/i} = \\left\[ \\begin{array}{ccc}\n  q\_{1}^{2}-q\_{2}^{2}-q\_{3}^{2}+q\_{4}^{2} & 2\\left(q\_{1}q\_{2}+q\_{3}q\_{4}\\right) & 2\\left(q\_{1}q\_{3}-q\_{2}q\_{4}\\right) \\\\\n2\\left(q\_{1}q\_{2}-q\_{3}q\_{4}\\right) & -q\_{1}^{2}+q\_{2}^{2}-q\_{3}^{2}+q\_{4}^{2} & 2\\left(q\_{2}q\_{3}+q\_{1}q\_{4}\\right)\\\\\n2\\left(q\_{1}q\_{3}+q\_{2}q\_{4}\\right) & 2\\left(q\_{2}q\_{3}-q\_{1}q\_{4}\\right) & -q\_{1}^{2}-q\_{2}^{2}+q\_{3}^{2}+q\_{4}^{2}\n\\end{array}\\right\]$
%[text]  なので，回転行列からクォータニオンを逆算する方法もあるが，少しtechnicalなことが必要．
%[text] cf. Markley, F. Landis. “Unit Quaternion from Rotation Matrix.” Journal of Guidance, Control, and Dynamics, vol. 31, no. 2, 2008, pp. 440-42.
%[text] cf. yoshimuLibrary/attitude/dcm2q.mlx,  yoshimuLibrary/attitude/dcm2qC.mlx
%[text] ## オイラー角の回転をそれぞれクォータニオンに対応させる方法
%[text] この方法だとオイラー角の時間履歴をfor文を使わずに一気に変換できる．
% q4がスカラー部の場合を想定
n = length(phi); % 3-2-1オイラー角の時間履歴データの計算にも対応させる．

q_phi = [zeros(n,2), sin(phi/2), cos(phi/2)]; % z軸まわり，nx4 matrix
q_theta = [zeros(n,1), sin(theta/2), zeros(n,1), cos(theta/2)];
q_psi = [sin(psi/2), zeros(n,2), cos(psi/2)];

%[text] クォータニオン積（quaternion multiplication）
%[text] $\\odot$と$\\otimes$の定義があるが，ここでは$\\otimes\n$を使用．（$\\odot$の場合は順番が変わる）
%[text] cf. [処方箋](https://ja.overleaf.com/project/592e7ff9ac6e3cde785ab97d)
q = qMult(4,1, q_psi, qMult(4,1, q_theta, q_phi)) %[output:526a50ea]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:526a50ea]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"q","rows":1,"type":"double","value":[["0.1436","0.1060","0.0343","0.9833"]]}}
%---
