%% propagateOrbit - 軌道要素から位置ベクトルを計算
%
% 入力:
%   oe: 軌道要素 [a, e, i, RAAN, omega, f], 1x6
%   tSpan: 時間配列 [s], 1xN
%   mu: 重力定数 [km^3/s^2]
%
% 出力:
%   rI: 慣性座標系での位置 [km], Nx3
%   vI: 慣性座標系での速度 [km/s], Nx3

function [rI, vI] = propagateOrbit(oe, tSpan, mu)

nPoints = length(tSpan);
n = sqrt(mu / oe(1)^3);

M = meanAnomaly(oe(2), oe(6)) + n * tSpan;
M = mod(M, 2*pi);

oeArray = repmat(oe, nPoints, 1);
oeArray(:, 6) = M';

[rI, vI] = oe2rv(oeArray, 0, mu);

end
