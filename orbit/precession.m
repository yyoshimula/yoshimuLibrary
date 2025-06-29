%[text] # precession from mean to TOD
%[text] `jd0`: Julian day, day, scalar
%[text] `jd1`: Julian day, day, scalar
%[text] `const: orbital constants`
%[text] zeta, z, theta: precession angles, rad
%[text] eta: angle of inclination between the two ecliptics.
%[text] Pi\_: angle from initial equinox to intersection of ecliptics,  rad
%[text] p: combined precession in longitude, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] \[1\] Montenbruck, O., and Gill, E., Satellite Orbits, Berlin, Heidelberg: Springer Science & Business Media, 2012. p.176
%[text] \[2\] Jean Meeus, "Astronomical Algorithms, 2nd Edition", pp.134-138
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe, nutationQ, precession.
function [zeta, z, theta, eta, Pi_, p] = precession(jd0, jd1, const)
% arguments
%     jd0 (:,1) {mustBeNumeric}
%     jd1 (:,1) {mustBeNumeric}
%     const
% end

T = (jd0 - const.J2000) ./ 36525.0;
t = (jd1 - jd0) ./ 36525.0;

zeta = (2306.2181 + 1.39656 .* T -0.000139 .* T.^2 ) .* t ...
    + (0.30188 - 0.000344 .* T) .* t.^2 ...
    + 0.017998 .* t.^3;
zeta = arcs2rad(zeta);

z = (2306.2181 + 1.39656 .* T - 0.000139 .* T.^2 ) .* t ...
    + (1.09468 + 0.000066 .* T ) .* t.^2 ...
    + 0.018203 .* t.^3;
z = arcs2rad(z);

theta = (2004.3109 - 0.85330 .* T - 0.000217 .* T.^2) .* t ...
    - (0.42665 + 0.000217 .* T) .* t.^2 ...
    - 0.041833 .* t.^3;
theta = arcs2rad(theta);

% angle of inclination between the two ecliptics. Eq. (21.5) [ref.2]
eta =(47.0029 - 0.06603 .* T + 0.000598 .* T.^2) .* t ...
    + (-0.03302 + 0.000598 .* T ) .* t.^2 ...    
    + 0.000060 .* t.^3;
eta = arcs2rad(eta);

% initial equinox to intersection of ecliptics, Eq. (21.5) [ref.2]
Pi_ = 3289.4789 .* T + 0.60622 .* T.^2 - (869.8089 + 0.50491 .* T) .* t ...    
    + 0.03536 .* t.^2;
Pi_ = deg2rad(174.876384) + arcs2rad(Pi_);

% combined precession in longitude, rad, Eq. (21.5) [ref.2]
p = (5029.0966 + 2.222226 .* T -0.0000042 .* T.^2) .* t ...
    + (1.11113 - 0.000042 .* T ) .* t.^2 - 0.000006 .* t.^3;
p = arcs2rad(p);
end

%[appendix]{"version":"1.0"}
%---
