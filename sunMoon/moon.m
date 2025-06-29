%[text] # Moon position w.r.t. J2000.0 frame.
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] const: orbital constants
%[text] ELP: ELP coefficients
%[text] ## outputs
%[text] `moonPos`: Moon position w.r.t. J2000.0 frame, km, nx3 matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] M. Chapront-Touzé and J. Chapront. The lunar ephemeris,  ELP 2000. Astronomy and Astrophysics, vol. 124, 1983, pp. 50-62.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sun.
function moonPos = moon(jd, const, ELP)
% arguments
%     jd (:,1) {mustBeNumeric}
%     const
%     ELP
% end

%[text] ## moon geocentric longitude, latitude, and distance w.r.t. mean ecliptic and equinox of J2000.0 frame
[lonM, latM, rM] = moonLonLatR(jd, const, ELP);
tmp = rM .* [cos(latM) * cos(lonM), cos(latM) * sin(lonM), sin(latM)];
%[text] ## conversion
%[text] ${\\bf r}\_{\\rm i, moon} = R\_{3}(-\\epsilon\_0) {\\bf r}\_{\\rm moon}\n$
% moon position vector at inertial frame (J2000.0 fram), km
tmpQ = [sin(-const.EPS0/2) 0 0 cos(-const.EPS0/2)];
tmpQ = repmat(tmpQ, size(tmp,1), 1);
moonPos = qRotation(4, tmp, tmpQ); % nx3

end


%[appendix]{"version":"1.0"}
%---
