%[text] # Sun's geocentric longitude, latitude, and distance for a given Julian date referred to the mean ecliptic and equinox of date
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] const: orbital constants
%[text] earthVSOP, earth VSOP87 constants
%[text] ## outputs
%[text] lon: Sun's geocentric longitude, rad
%[text] lat: Sun's geocentric latitude, rad
%[text] r: Sun's geocentric distance, AU
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Jean Meeus, "Astronomical Algorithms, 2nd edition", p.166 and p. 217.  
%[text] To obtain the geocentric longitude and latitude of the Sun, add 180 (deg) to Earth's heliocentric longitude and change the sign of Earth's heliocentric latitude.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, vsopConst, precession.
function [lon, lat, r] = sunLonLatR(jd, const, earthVSOP)
% arguments
%     jd (:,1) {mustBeNumeric}
%     const
%     earthVSOP
% end

% earth's heliocentric longitude, latitude, and distance
[lon, lat, r] = earthVSOP87(jd, earthVSOP);

lon = lon + pi;
lat = -1 .* lat;

[~, ~, ~, eta, Pi_, p] = precession(const.J2000, jd, const);

A = sin(eta) .* sin(lat) + cos(eta) .* cos(lat) .* sin(p + Pi_ - lon);
B = cos(lat) .* cos(p + Pi_ - lon);
C = cos(eta) .* sin(lat) - sin(eta) .* cos(lat) .* sin(p + Pi_ - lon);

lon = Pi_ - atan2(A,B);
lat = asin(C);

end

%[appendix]{"version":"1.0"}
%---
