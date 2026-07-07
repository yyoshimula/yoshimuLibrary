%[text] # Compute derived orbital state quantities from orbital elements
%[text] 軌道要素から派生量（p, h, n, u, M, uM, r）と位置・速度ベクトルを計算
%[text] ## inputs
%[text] `oe`: struct with fields a, e, inc, raan, w, nu (semi-major axis, eccentricity, inclination, RAAN, argument of perigee, true anomaly)
%[text] `mu`: gravitational constant (unit consistent with oe.a)
%[text] ## outputs
%[text] `oe`: input struct augmented with the following fields
%[text] `oe.p`: semilatus rectum
%[text] `oe.h`: orbital angular momentum
%[text] `oe.n`: mean motion
%[text] `oe.u`: true argument of latitude (mod 2*pi)
%[text] `oe.M`: mean anomaly (mod 2*pi)
%[text] `oe.uM`: mean argument of latitude (mod 2*pi)
%[text] `oe.r`: radius
%[text] `oe.rVec`: position vector
%[text] `oe.vVec`: velocity vector
%[text] See also oe2rv, meanAnomaly, gve.
function oe = calcOrbitalState(oe, mu)
% Calculate orbital parameters
oe.p = oe.a * (1 - oe.e^2); % semilatus rectum
oe.h = sqrt(mu * oe.p); % orbital angular momentum
oe.n = sqrt(mu / oe.a^3); % mean motion
oe.u = oe.w + oe.nu; % true argument of latitude
oe.M = meanAnomaly(oe.e, oe.nu);
oe.uM = oe.w + oe.M; % mean argument of latitude
oe.r = oe.p / (1.0 + oe.e * cos(oe.nu)); % radius

oe.u = mod(oe.u, 2*pi);
oe.M = mod(oe.M, 2*pi);
oe.uM = mod(oe.uM, 2*pi);

% Calculate position and velocity vectors
[rTmp, vTmp] = oe2rv([oe.a, oe.e, oe.inc, oe.raan, oe.w, oe.nu], 1, mu);

oe.rVec = rTmp(:);
oe.vVec = vTmp(:);

end

%[appendix]{"version":"1.0"}
%---
