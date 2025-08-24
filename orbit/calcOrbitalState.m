
function oe = calcOrbitalState(oe, mu)
% Calculate orbital parameters
oe.p = oe.a * (1 - oe.e^2); % semilatus rectum
oe.h = sqrt(mu * oe.p); % orbital angular momentum
oe.n = sqrt(mu / oe.a^3); % mean motion
oe.u = oe.ome + oe.f; % true argument of latitude
oe.M = meanAnomaly(oe.e, oe.f);
oe.uM = oe.ome + oe.M; % mean argument of latitude
oe.r = oe.p / (1.0 + oe.e * cos(oe.f)); % radius

% Calculate position and velocity vectors
[rTmp, vTmp] = oe2rv([oe.a, oe.e, oe.inc, oe.raan, oe.ome, oe.f], 1, mu);

oe.rVec = rTmp(:);
oe.vVec = vTmp(:);

end
