%[text] # Converting TOD frame to MOD frame
%[text] True Equator Mean EquinoxをMean equinox of date座標へ変換（IAU-1980 theory of nutation.）
%[text] `jd`: julian day, day
%[text] `i`: inclination, rad
%[text] `Ome`: longitude of the ascending node, rad
%[text] `w`: argument of the perihelion, rad
%[text] `iMod`: inclination at MOD, rad
%[text] `OmeMod`: longitude of the ascending node at MOD, rad
%[text] `wMod`: argument of the perihelion at MOD, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 3rd ed.," pp.229-231.
%[text] Jean Meeus, "Astronomical Algorithms, 2nd ed.," pp.143-148. 
%[text] ## revisions
%[text] 20210601  y.yoshimura
%[text] See also teme2Mod, mod2J2000, obliquity, nutation.
function [iMod, OmeMod, wMod] = tod2Mod(jd, i, Ome, w, const)

%[text] $\\Delta \\Psi, \\Delta\\epsilon$
e0 = obliquity(jd);
DCM = nutationDCM(jd, const);

tmp = [sin(Ome) * sin(i)
    -cos(Ome) * sin(i)
    cos(i)];

iMod = acos(DCM(:,3)' * tmp);

OmeMod = atan2(DCM(:,1)' * tmp, DCM(:,2)' * (-tmp));

wMod = w - atan2(-cos(Ome) * DCM(1,3) - sin(Ome) * DCM(2,3), ...
    -sin(Ome) * cos(i) * DCM(1,3) + cos(Ome) * cos(i) * DCM(2,3) + sin(i) * DCM(3,3));

%[text] $\\Omega, \\omega \\in \[0, 2\\pi\]$
wMod = mod(wMod, 2*pi);
OmeMod = mod(OmeMod, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
