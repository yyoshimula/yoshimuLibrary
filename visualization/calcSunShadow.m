%% calcSunShadow - 太陽位置と日食判定
%
% 入力:
%   jdArray: Julian day配列, Nx1 or 1xN
%   rI: 慣性座標系での衛星位置 [km], Nx3
%   const: 軌道定数
%
% 出力:
%   sunPos: 太陽位置 [km], Nx3
%   nu: 日照係数 (1:日照, 0:影), Nx1

function [sunPos, nu] = calcSunShadow(jdArray, rI, const)

earthVSOP = vsopConst();
sunPos = sun(jdArray(:), const, earthVSOP);
nu = shadow(rI, sunPos, const.RS, const.RE);

end
