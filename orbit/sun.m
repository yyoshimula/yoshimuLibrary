%[text] # Sun position w.r.t. J2000.0 frame.
%[text] ## inputs
%[text] `jd`: Julian date of UTC, day, nx1 vector
%[text] const: orbital constants
%[text] earthVSOP: earth VSOP coefficients
%[text] ## outputs
%[text] `sunPos`: Sun position w.r.t. J2000.0 frame, km, nx3 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sunLonLatR, earthVSOP, vsopConst.
function sunPos = sun(jd, const, earthVSOP)
% arguments
%     jd (:,1) {mustBeNumeric}
%     const
%     earthVSOP
% end

%[text] ## Sun geocentric longitude, latitude, and distance w.r.t. mean ecliptic and equinox of J2000.0 frame
[lonS, latS, sunAU] = sunLonLatR(jd, const, earthVSOP);
tmp = au2km(sunAU, const) .* [cos(latS) .* cos(lonS), cos(latS) .* sin(lonS), sin(latS)];
%[text] ## conversion
%[text] ${\\bf r}\_{\\rm i, sun} = R\_{3}(-\\epsilon\_0) {\\bf r}\_{\\rm sun}\n$
% moon position vector at inertial frame (J2000.0 fram), km
tmpQ = [sin(-const.EPS0/2) 0 0 cos(-const.EPS0/2)];
tmpQ = repmat(tmpQ, size(tmp,1), 1);
sunPos = qRotation(4, tmp, tmpQ); % nx3

end

%[appendix]{"version":"1.0"}
%---
