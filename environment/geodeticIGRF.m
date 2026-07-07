function b = geodeticIGRF(jd, lat, lon, alt, coefs)
% ----------------------------------------------------------------------
%   calculate geocentric magnetic vector with IGRF model
%    20210428  y.yoshimura
%    Inputs: jd: Julian day, day
%            lat: geodetic latitude, rad
%            lon: geodetic longitude, rad
%            alt: geodetic altitude, km
%   Outputs: b: magnetic field, [northward, eastward, downward], nT
%   related function files:
%   note:
%   cf:
%   revisions;
%   
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

[year, ~, ~, ~, ~, ~] = jd2gc(jd);
% decimal year for IGRF coefficient interpolation
decYear = year + (jd - gc2jd(year, 1, 1, 0, 0, 0)) ...
                 ./ (gc2jd(year + 1, 1, 1, 0, 0, 0) - gc2jd(year, 1, 1, 0, 0, 0));

% NOTE (20260707): 旧実装は外部関数 igrfs.m（リポジトリ外・未同梱）に依存しており
%   実行時に "Undefined function 'igrfs'" で失敗していた。リポジトリ内蔵の igrf12
%   （igrfsyn12 mex + igrf12coeffs.txt）へ置き換え。lat/lon は rad のまま渡す
%   （igrf12 が内部で deg 変換する）。引数 coefs は igrf12 が独自の IGRF-12 係数を
%   用いるため使用しない（後方互換のため引数だけ残置）。
b = igrf12(decYear, alt, lat, lon); % 1x3, nT, NED (northward, eastward, downward)
end