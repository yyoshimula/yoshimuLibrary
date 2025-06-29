%[text] # precession direction cosine matrix (DCM) from J2000(ICRF) to mean date
%[text] `jd0`: Julian day, day, scalar
%[text] `jd1`: Julian day, day, scalar
%[text] const: orbital constants
%[text] dcm: direction cosine matrix from jd0 to jd1 when jd0 = J2000.0 means the rotation from J2000(ICRF) to MOD frame
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Montenbruck, O., and Gill, E., Satellite Orbits, Berlin, Heidelberg: Springer Science & Business Media, 2012. p.176
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also oe2roe, nutationQ, precession.
function dcm = precessionDCM(jd0, jd1, const)
% arguments
%     jd0 (:,1) {mustBeNumeric}
%     jd1 (:,1) {mustBeNumeric}    
%     const
% end

[zeta, z, theta, ~, ~, ~] = precession(jd0, jd1, const);

% dcm = dcm1axis(3, -pi/2 - z) * dcm1axis(1, theta) * dcm1axis(3, pi/2 - zeta);
dcm = dcm1axis(3, -z) * dcm1axis(2, theta) * dcm1axis(3, -zeta);

end

%[appendix]{"version":"1.0"}
%---
