%[text] # Converting TEME frame to J2000 frame
%[text] True Equator Mean Equinox (TEME)をJ2000座標（FK5）へNewcomb's precessional constantsを用いて変換
%[text] ## inputs
%[text] `jd`: julian day, day
%[text] `i`: inclination, rad
%[text] `Ome`: longitude of the ascending node, rad
%[text] `w`: argument of the perihelion, rad
%[text] ## outputs
%[text] `iJ`: inclination at J2000.0, rad
%[text] `wJ`: argument of the perihelion at J2000.0, rad
%[text] `OmeJ`: longitude of the ascending node at J2000.0, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 3rd ed.," pp.229-231.
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th ed.," pp.231-231.
%[text] Jean Meeus, "Astronomical Algorithms, 2nd ed.," pp.143-148. 
%[text] ## revisions
%[text] 2021020209  y.yoshimura
%[text] See also teme2Mod, mod2J2000, mean2Osc.
function [iJ, OmeJ, wJ] = teme2J2000(jd, i, Ome, w, const)

%[text] ### TEME to TOD
epsi = obliquity(jd);
%[text] $\\Delta \\Psi, \\Delta\\epsilon$
[dPsi, dEpsi] = nutation(jd, const);
Ome = Ome + dPsi .* cos(epsi + dEpsi);

%[text] ### TOD → MOD → J2000
[iMod, OmeMod, wMod] = tod2Mod(jd, i, Ome, w, const);
[iJ, OmeJ, wJ] = mod2J2000(jd, iMod, OmeMod, wMod, const);

end

%[appendix]{"version":"1.0"}
%---
