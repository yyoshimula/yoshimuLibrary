%[text] # calculating solar radiation pressure (SRP) forces using cannonball model
%[text] ## inputs
%[text] `satAm`: area-to-mass ratio
%[text] `sunRel:` sun vector from satellite to Sun expressed with inertial frame, nx3 unit vector
%[text] `d`: distance between satellite and Sun, m
%[text] `const`: orbital constants
%[text] Cr: (optional), reflectivity, default value is 1.2
%[text] ## output
%[text] `aSRP`: SRP acceleration w.r.t. inertial frame, $\\rm m/s^2$, nx3 matrix
%[text] ## note
%[text] NA
%[text] ## references
%[text] Montenbruck, Oliver, & Eberhard Gill. Satellite Orbits. Springer Science & Business Media, 2012., p79
%[text] ## revisions
%[text] 20221010  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readSC, orbitConst.
function aSRP = srpCannon(satAm, sunRel, d, const, Cr)
% arguments (Input)
%     satAm
%     sunRel (:,3) {mustBeNumeric}
%     d (:,1) {mustBeNumeric}
%     const
%     Cr (1,1) {mustBeNumeric} = 1.2
% end
% arguments (Output)
%     aSRP (:,3)
% end

if nargin < 5
    Cr = 1.2;
end

%[text] ## parameters
dAU = km2au(d ./ 10^3, const); % AU

coeff = -const.S0 / const.c * Cr * satAm;

%[text] ## SRP
sunRel = sunRel ./ vecnorm(sunRel, 2, 2);

aSRP = coeff ./ dAU.^2 .* sunRel; % nx3 matrix

end

%[appendix]{"version":"1.0"}
%---
