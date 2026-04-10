function [F10, F10a, Kp] = lookupSolarGeoIndex(jd, sw)
% LOOKUPSOLARGEOINDEX  CelesTrak SW-All.csv から太陽・地磁気指数を取得 (JR1971用)
%
% Inputs:
%   jd  - Julian day (UTC)
%   sw  - 以下のいずれか:
%         (a) loadSpaceWeather() で読み込んだ構造体 (推奨)
%         (b) SW-All.csv のファイルパス（文字列）。初回読み込み後は persistent
%             でキャッシュされ、同一パスなら再読み込みされない。
%
% Outputs:
%   F10  - 10.7 cm solar flux (1日ラグ, 観測値) [sfu]
%   F10a - 81-day centered average F10.7 (観測値) [sfu]
%   Kp   - Kp 指数 (jd 時刻の 3時間ブロック値, ラグなし)
%
% Data source:
%   CelesTrak Space Weather file SW-All.csv
%   https://celestrak.org/SpaceData/SW-All.csv
%
% Dependencies:
%   yoshimuLibrary: jd2gc, gc2jd
%
% Example:
%   [F10, F10a, Kp] = lookupSolarGeoIndex(gc2jd(2024,1,1,12,0,0), 'SW-All.csv');

%% 入力解決：パス文字列ならキャッシュから読み込む
persistent cachedSW cachedPath
if ischar(sw) || isstring(sw)
    swPath = char(sw);
    if isempty(cachedSW) || ~strcmp(cachedPath, swPath)
        cachedSW   = loadSpaceWeather(swPath);
        cachedPath = swPath;
    end
    sw = cachedSW;
elseif ~isstruct(sw)
    error('lookupSolarGeoIndex:BadInput', ...
        'sw must be a struct from loadSpaceWeather or a path to SW-All.csv.');
end

%% F10 / F10a (1日ラグ, 観測値)
% floor(jd - 0.5) は jd の UTC 暦日における「前日の正午 JD (整数)」を返す。
JD_lag = floor(jd - 0.5);
iLag = round(JD_lag - sw.jdNoonStart) + 1;
if iLag < 1 || iLag > sw.nDays
    error('lookupSolarGeoIndex:OutOfRange', ...
        'JD %.3f (lag date %s) is outside SW-All range [%s, %s].', ...
        jd, jdToDateStr(JD_lag), ...
        jdToDateStr(sw.jdNoonStart), jdToDateStr(sw.jdNoonEnd));
end
F10  = sw.F107obs(iLag);
F10a = sw.F107obsCenter81(iLag);

if isnan(F10) || isnan(F10a)
    error('lookupSolarGeoIndex:MissingData', ...
        'F10/F10a is NaN at %s (likely beyond observed/centered-81 coverage).', ...
        jdToDateStr(JD_lag));
end

%% Kp (jd 時刻の 3時間ブロック, ラグなし)
% floor(jd + 0.5) は jd の UTC 暦日における「当日の正午 JD (整数)」を返す。
JD_today = floor(jd + 0.5);
iNow = round(JD_today - sw.jdNoonStart) + 1;
if iNow < 1 || iNow > sw.nDays
    error('lookupSolarGeoIndex:OutOfRange', ...
        'JD %.3f is outside SW-All range [%s, %s].', ...
        jd, jdToDateStr(sw.jdNoonStart), jdToDateStr(sw.jdNoonEnd));
end
[~, ~, ~, hr, ~, ~] = jd2gc(jd);
block = floor(hr / 3) + 1;             % 1..8
Kp = sw.Kp(iNow, block);
end

%% =========================================================================
%  Local Functions
%  =========================================================================

function s = jdToDateStr(jdNoon)
% 整数正午 JD を 'yyyy-mm-dd' 形式に整形
    [y, m, d] = jd2gc(jdNoon);
    s = sprintf('%04d-%02d-%02d', y, m, d);
end
