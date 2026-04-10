function sw = loadSpaceWeather(csvPath)
% LOADSPACEWEATHER  CelesTrak SW-All.csv を読み込み構造体を返す
%
% Input:
%   csvPath - SW-All.csv のパス
%
% Output:
%   sw.jdNoonStart     - 先頭行の正午 JD (整数)
%   sw.jdNoonEnd       - 末尾行の正午 JD (整数)
%   sw.nDays           - 行数 (= 連続日数)
%   sw.Kp              - Kp 3時間ブロック (N x 8), 実数
%   sw.Ap              - Ap 3時間ブロック (N x 8)
%   sw.ApAvg           - 日平均 Ap (N x 1)
%   sw.F107obs         - 観測 F10.7 (N x 1) [sfu]
%   sw.F107obsCenter81 - 81日中央平均 F10.7 観測値 (N x 1) [sfu]
%
% Data source:
%   CelesTrak Space Weather file SW-All.csv
%   https://celestrak.org/SpaceData/SW-All.csv
%
% Dependencies:
%   yoshimuLibrary: gc2jd

if exist(csvPath, 'file') ~= 2
    error('loadSpaceWeather:FileNotFound', ...
        'Space weather file not found: %s', csvPath);
end

% --- DATE 列のみテキストとして読み込む ------------------------------------
opts = detectImportOptions(csvPath, 'Delimiter', ',');
if isprop(opts, 'VariableNamingRule')
    opts.VariableNamingRule = 'preserve';
end
opts = setvartype(opts, 'DATE', 'char');
T = readtable(csvPath, opts);

sw.nDays = height(T);

% 先頭・末尾の日付を yyyy-mm-dd から (Y, M, D) に分解 → gc2jd で正午 JD
[y0, m0, d0] = parseYmd(T.DATE{1});
[yN, mN, dN] = parseYmd(T.DATE{end});
sw.jdNoonStart = gc2jd(y0, m0, d0, 12, 0, 0);
sw.jdNoonEnd   = gc2jd(yN, mN, dN, 12, 0, 0);

% 連続日数チェック
expectedDays = round(sw.jdNoonEnd - sw.jdNoonStart) + 1;
if expectedDays ~= sw.nDays
    warning('loadSpaceWeather:NonConsecutive', ...
        'SW-All.csv の行数 (%d) と日付スパン (%d) が一致しません。', ...
        sw.nDays, expectedDays);
end

% --- Kp / Ap / F10.7 列の取り出し ----------------------------------------
% CelesTrak Kp は Kp*10 で格納 → 10で割って実数化
sw.Kp = double([T.KP1, T.KP2, T.KP3, T.KP4, T.KP5, T.KP6, T.KP7, T.KP8]) / 10;
sw.Ap = double([T.AP1, T.AP2, T.AP3, T.AP4, T.AP5, T.AP6, T.AP7, T.AP8]);
sw.ApAvg = double(T.AP_AVG);

% 列名にドットを含むため curly-brace 形式で取得
sw.F107obs         = double(T{:, 'F10.7_OBS'});
sw.F107obsCenter81 = double(T{:, 'F10.7_OBS_CENTER81'});
end

function [y, m, d] = parseYmd(s)
% 'yyyy-mm-dd' 形式の文字列から年月日を抽出
    p = sscanf(s, '%d-%d-%d');
    y = p(1); m = p(2); d = p(3);
end
