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
