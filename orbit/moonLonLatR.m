%[text] # Moon's geocentric longitude, latitude, and distance w.r.t. the mean ecliptic and equinox of J2000.0.
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] const: orbital constants
%[text] ELP: ELP coefficients
%[text] ## outputs
%[text] `lonM`: Moon's geocentric longitude, rad
%[text] `latM`: Moon's geocentric latitude, rad
%[text] `rM`: Moon's geocentric distance, km
%[text] ## note
%[text] NA
%[text] ## references 
%[text] M. Chapront-Touzé and J. Chapront. The lunar ephemeris,  ELP 2000. Astronomy and Astrophysics, vol. 124, 1983, pp. 50-62.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sun.
function [lonM, latM, rM] = moonLonLatR(jd, const, ELP)
% arguments
%     jd (:,1) {mustBeNumeric}
%     const        
%     ELP
% end

%[text] ## moon position w.r.t. ELP 2000 frame
  [lon, lat, rM] = moonELP(jd, ELP);

  [~, ~, ~, eta, Pi_, p] = precession(const.J2000, jd, const);

%[text] ## conversion
  A = sin(eta) * sin(lat) + cos(eta) * cos(lat) * sin(p + Pi_ - lon);
  B = cos(lat) * cos(p + Pi_ - lon);
  C = cos(eta) * sin(lat) - sin(eta) * cos(lat) * sin(p + Pi_ - lon);
  
  lonM = Pi_ - atan2(A, B);  
  latM = asin(C);  

end

%[appendix]{"version":"1.0"}
%---
