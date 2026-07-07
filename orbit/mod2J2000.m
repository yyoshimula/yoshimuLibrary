%[text] # Converting MOD frame to J2000 frame
%[text] Mean equinox of date (MOD)座標をJ2000座標へ変換（using Newcomb's precessional constants）
%[text] ## inputs
%[text] `jd`: julian day, day
%[text] `i`: inclination, rad
%[text] `Ome`: longitude of the ascending node, rad
%[text] `w`: argument of the perihelion, rad
%[text] `const`: constant parameters for orbit propagation
%[text] ## outputs
%[text] `iJ`: inclination at J2000.0, rad
%[text] `wJ`: argument of the perihelion at J2000.0, rad
%[text] `OmeJ`: longitude of the ascending node at J2000.0, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 3rd ed.," pp.229-231.
%[text] Jean Meeus, "Astronomical Algorithms, 2nd ed.," pp.143-148. 
%[text] ## revisions
%[text] 20210601  y.yoshimura
%[text] See also teme2J2000, orbitConst.
function [iJ, OmeJ, wJ] = mod2J2000(jd, i, Ome, w, const)
% arguments
%     jd (:,1) {mustBeNumeric}
%     i (:,1) {mustBeNumeric}
%     Ome (:,1) {mustBeNumeric}
%     w (:,1) {mustBeNumeric}
%     const
% end

[zeta, z, theta, ~, ~, ~] = precession(const.J2000, jd, const);

iJ = acos(cos(theta) * cos(i) - sin(theta) * sin(Ome - z) * sin(i));

OmeJ = atan2(cos(theta) * sin(Ome - z) * sin(i) ...
    + sin(theta) * cos(i), cos(Ome - z) * sin(i)) - zeta;

wJ = w - atan2(sin(theta) * cos(Ome - z), ...
    cos(theta) * sin(i) + sin(theta) * sin(Ome - z) * cos(i));

%[text] $\\Omega, \\omega \\in \[0, 2\\pi\]$
OmeJ = mod(OmeJ, 2*pi);
wJ = mod(wJ, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
