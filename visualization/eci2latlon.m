%% eci2latlon - ECI座標から緯度経度を計算
%
% 入力:
%   rI: 慣性座標系での位置 [km], Nx3
%   GMSTarray: Greenwich Mean Sidereal Time [rad], Nx1
%   const: 軌道定数
%
% 出力:
%   lat: 測地緯度 [rad], Nx1
%   lon: 経度 [rad], Nx1

function [lat, lon] = eci2latlon(rI, GMSTarray, const)

nPoints = size(rI, 1);
rECEF = zeros(size(rI));

for k = 1:nPoints
    rECEF(k,:) = (dcm1axis(3, -GMSTarray(k)) * rI(k,:)')';
end

[lat, lon, ~] = ecef2LatLonH(rECEF, const);

end
