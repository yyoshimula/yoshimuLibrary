%[text] # Associated Legendre多項式
%[text] P_l^m(x) を再帰関係式により計算する（Condon-Shortley 位相込み・非正規化）
%[text] ### 非正規化Associated Legendre多項式（全位数）
function [P, m_values] = associatedLegendre(l, x) %#codegen
    % 全ての位数 m = 0, 1, ..., l に対してP_l^m(x)を計算
    % 入力:
    %   l: 次数 (非負整数)
    %   x: 評価点 (スカラーまたはベクトル)
    % 出力:
    %   P: (l+1) × length(x) の行列、P(i,j) = P_l^{m_i-1}(x_j)
    %   m_values: 位数のベクトル [0, 1, ..., l]
    
    m_values = 0:l;
    P = zeros(length(m_values), length(x));
    
    % 方法1: 再帰関係式を使用
    for idx = 1:length(m_values)
        m = m_values(idx);
        P(idx, :) = legendreRecursive(l, m, x);
    end
end
%[text] ## 再帰関係式による実装
function P = legendreRecursive(l, m, x)
    % Associated Legendre多項式 P_l^m(x) を再帰的に計算
    % 入力:
    %   l: 次数 (非負整数)
    %   m: 位数 (0 <= m <= l)
    %   x: 評価点 (スカラーまたはベクトル)
    
    % 入力チェック
    % if m < 0 || m > l
    %     error('mは0以上l以下である必要があります');
    % end
    
    % if any(abs(x) > 1)
    %     warning('|x| > 1 の点が含まれています');
    % end

    P = 0; % initialize P
    
    % P_m^m(x) = (-1)^m * (2m-1)!! * (1-x^2)^(m/2)
    if l == m
        P = (-1)^m * double_factorial(2*m-1) * (1-x.^2).^(m/2);
        return;
    end
    
    % P_{m+1}^m(x) = x * (2m+1) * P_m^m(x)
    if l == m + 1
        Pmm = (-1)^m * double_factorial(2*m-1) * (1-x.^2).^(m/2);
        P = x .* (2*m + 1) .* Pmm;
        return;
    end
    
    % 再帰関係式を使用
    % (l-m)P_l^m = x(2l-1)P_{l-1}^m - (l+m-1)P_{l-2}^m
    Pmm = (-1)^m * double_factorial(2*m-1) * (1-x.^2).^(m/2);
    Pmp1m = x .* (2*m + 1) .* Pmm;
    
    for k = m+2:l
        P = (x .* (2*k-1) .* Pmp1m - (k+m-1) .* Pmm) / (k-m);
        Pmm = Pmp1m;
        Pmp1m = P;
    end
end
%%
%[text] ## 補助関数: 二重階乗
function result = double_factorial(n)
    % n!! を計算
    if n <= 0
        result = 1;
    elseif mod(n, 2) == 0
        % 偶数の場合: n!! = n * (n-2) * (n-4) * ... * 2
        result = prod(2:2:n);
    else
        % 奇数の場合: n!! = n * (n-2) * (n-4) * ... * 1
        result = prod(1:2:n);
    end
end

%[appendix]{"version":"1.0"}
%---
