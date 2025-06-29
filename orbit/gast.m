%[text] # calculating Greenwich apparent sidereal time (GAST) from Julian day
%[text] calculation based on IAU-76/FK5
%[text] ## inputs
%[text] `jd`: Julian day, day, nx1 vector
%[text] `const`: orbital constants
%[text] ## output
%[text] `GAST`: Greenwich apparent sidereal time, rad, nx1 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Jean Meeus, "Astronomical Algorithms, 2nd edition", p.87. Eq.(12.3) or (12.4)
%[text] ## revisions
%[text] 20160419  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gmst.
function GAST = gast(jd, const)
% arguments
%     jd (:,1) {mustBeNumeric}    
%     const
% end

GMST = gmst(jd); % Greenwich mean sidreal time
meanEpsi = obliquity(jd); % mean obliquity of the ecliptic, rad

[dPsi, dEpsi] = nutation(jd, const);

GAST = GMST + dPsi .* cos(meanEpsi + dEpsi);

end

%[appendix]{"version":"1.0"}
%---
