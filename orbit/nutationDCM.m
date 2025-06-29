%[text] # nutation rotation matrix from mean of date (MOD) to true of date(TOD)
%[text] # (IAU-76/FK5)
%[text] ## inputs
%[text] `jd`: Julian day, day, scalar
%[text] `const:` orbital constants
%[text] ## output
%[text] dcm: direction cosine matrix from MOD to TOD
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Explanatory Supplement to the Astronomical Almanac, p. 114. 
%[text] Satellite Orbits, p.180
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe, orbitConst.
function dcm = nutationDCM(jd, const)
% arguments
%     jd (:,1) {mustBeNumeric}
%     const
% end

%[text] ## mean obliquity of the ecliptic at the epoch, rad
e = obliquity(jd);

% dPsi: nutation in longitude, rad
% dEpsi: nutation in obliquity, rad
[dPsi, dEpsi] = nutation(jd, const);

ePrime = dEpsi + e;
%[text] R1 = DCM1axis(1, -e\_prime);
%[text] R2 = DCM1axis(3, -dPsi);
%[text] R3 = DCM1axis(1, e);
%[text] dcm = R1 \* R2 \* R3;
%[text] ## from mean of date to true of date
dcm = dcm1axis(1, -ePrime) * dcm1axis(3, -dPsi) * dcm1axis(1, e);


end

%[appendix]{"version":"1.0"}
%---
