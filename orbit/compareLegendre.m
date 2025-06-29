%[text] # comparing associated Legendre function calculations
%[text] このcodingだと，degがとても大きいとき（150以上くらい？）以外は 1.の書き方が早い
clc
clear
cls
%[text] ## condition setting
deg = 10;

phi = deg2rad(30);
sp = sin(phi);

P = zeros(deg+1, deg+1); % matlabはindexが1始まりなので+1

%[text] ## 1. MATLABのlegendre関数を使う方法
tic 
for i = 1:deg
    tmp = 0:i; % malabのlegendre functionは(-1)^mが付くため
    P(i+1,1:i+1) = (-1).^tmp .* legendre(i, sp)';
end

toc %[output:3f976e52]
P
%[text] ## 2.Recursiveに計算する方法
P = zeros(deg+1, deg+1); % matlabはindexが1始まりなので+1
tic

P(1,1) = 1;
P(2,1) = (2 * 0 + 1) * sp * P(1,1);
P(2,2) = (2 * 1 - 1) * sqrt(1 - sp^2) * P(1,1);
for n = 2:deg    
    P(n+1,n+1) = (2 * n - 1) * sqrt(1 - sp^2) * P(n,n);
    P(n+2,n+1) = (2 * n + 1) * sp * P(n+1,n+1);
    for m = 0:n-1
        P(n+1,m+1) =  1 / (n - m) * ((2*n - 1) * sp * P(n,m+1) - (n + m - 1) * P(n-1,m+1));        
    end
end
toc %[output:73d8962a]
P




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":19.4}
%---
%[output:3f976e52]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.001105 秒です。\n","truncated":false}}
%---
%[output:73d8962a]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.003852 秒です。\n","truncated":false}}
%---
