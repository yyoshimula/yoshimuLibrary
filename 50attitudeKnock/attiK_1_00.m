%[text] ## 姿勢計算50本ノック
%[text] ## 1-00
%[text] 適当な慣性テンソルJを設定して，その主慣性モーメントを求めよ
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
%[text] 固有ベクトル`vと`固有値`d`を求める`eig`関数を使うのが簡単
%[text] 主慣性モーメント = 慣性テンソルの固有値なので，dがそのまま主慣性モーメントになる．．
%[text] ただし，三角不等式を満たさないといけない．つまり次式全てを満たさないと物理的におかしい慣性テンソルになっている可能性がある.
%[text] 主慣性モーメントを$J\_x,J\_y,J\_z\n$とすると
%[text] $J\_x+J\_y \> J\_z\\\\\nJ\_x + J\_z \> J\_y \\\\\nJ\_y + J\_z \> J\_x\n$
%[text] を満たす．
[v, d] = eig(J);
d

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
