% Converting absolute orbital elements to position and velocity vector
% inputs
% oe: orbital elements  (m or km, -, rad, rad, rad, rad)
% flag: 1:= true anomaly, 0:= mean anomaly 
% GE: gravitational constant of the Earth, m or km (unit must be unified with the position, velocity, and semi-major axis)
% outputs
% r: position vector, nx3 vector, m or km 
% v: velocity vector, nx3 vector, m/s or km/s 
% note
% 重力定数と位置・速度，軌道長半径は単位を合わせること
% references 
% Vallado, D. A., and McClain, W. D., Fundamentals of Astrodynamics and Applications,revisions
% revisions
% 20231204 y.yoshimura, GEを直接引数として設定
% 20211027  y.yoshimura, y.yoshimula@gmail.com
% See also roe2losApprox, trueAnomaly.

function [r, v] = oe2rv_(oe, flag, GE)

%% orbital elements and constants
small = 1.0e-10;
mu = GE; % 地心重力定数 [km^3/s^2 or m^3/s^2]

a = oe(:,1);
e = oe(:,2);
inc = oe(:,3);
raan = oe(:,4);
ome = oe(:,5);
if flag == 1 % true anomaly
    f = oe(:,6);
else % mean anomaly
    [f, ~] = trueAnomaly_(a, e, mod(oe(:,6),2*pi));
end

% if orbit is circular, ome = 0.0;
ome = ome .* (e > small);

% if orbit is equatorial, raan = 0.0;
raan = raan .* (inc > small) .* (abs(inc - pi) > small);

p = a .* (1 - e.^2); % semi-latus rectum
p = 0.0001 .* (abs(p) < small) + p .* ( abs(p) >= small);

temp = p ./ (1.0  + e .* cos(f));

% position in orbital plane
rPQW = [temp .* cos(f), temp .* sin(f), zeros(length(f),1)];
vPQW = [-sin(f).*sqrt(mu)./sqrt(p), (e+cos(f)) .* sqrt(mu) ./ sqrt(p), zeros(length(f),1)];

% transformation to IJK by quaternion
q = zxz2q_(4, raan, inc, ome); % quaternion from IJK to PQW
q = qInv_(4, q);% quaternion from PQW to IJK

r = qRotation_(4, rPQW, q); % position@inertial frame, Nx3 matrix
v = qRotation_(4, vPQW, q); % velocity@inertial frame, Nx3 matrix

end
