%[text] ## 姿勢計算50本ノック
%[text] ## 1-05 適当な1000セットのクォータニオン（つまり，1000x4の行列）を準備し，${\\bf q} = \\left\[\\frac{1}{\\sqrt{4}}, \\frac{1}{\\sqrt{4}}, \\frac{1}{\\sqrt{4}}, \\frac{1}{\\sqrt{4}}\\right\]^T$とのクォータニオン積を計算せよ．
clc
clear
close all

q = ones(1,4) ./ sqrt(4); % column vectorで定義してもよい
%%
%[text] ## 適当な1000セットのクォータニオン
N = 1000;
p = rand(N, 4);
p = p ./ vecnorm(p, 2, 2); % vecnormで2乗ノルムを計算して正規化するのがポイント
pv = p(:,1:3); % pのvector part
%%
%[text] ## quaternion multiplication
%[text] クォータニオンは$\\otimes\n$と$\\odot$の2通りの定義がある．
%[text] 自身が使ってる定義の意味等を理解していればどっちでもok
%[text] `.*` や dot関数の次元（`dim`）をちゃんと指定することでfor文を使わなくて済む（地味に大事）
%[text] ### ${\\bf p} \\odot {\\bf q}$
Q = repmat(q, N,1); % qを1000個並べたmatrixを準備．
Qv = Q(:,1:3); % vector partが1から3と仮定

ans = [p(:,4) .* Qv + Q(:,4) .* pv + cross(pv, Qv), Q(:,4) .* p(:,4) - dot(pv, Qv, 2)]

%[text] ノルムが1か確認
vecnorm(ans, 2, 2)
%[text] ${\\bf p} \\otimes {\\bf q}$
ans = [p(:,4) .* Qv + Q(:,4) .* pv - cross(pv, Qv), Q(:,4) .* p(:,4) - dot(pv, Qv, 2)]

%[text] ノルムが1か確認
vecnorm(ans, 2, 2)
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
