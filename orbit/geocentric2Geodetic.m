%[text] # Converting geocentric rectangular coordinates to geodectic latitude, longitude
%[text] x, y, z; geocentric coordinates
%[text] a: Equatorial radius of reference ellipsoid. 
%[text] f: Flattening of reference ellipsoid: f = ( a - b ) / a, where a and b are ellipsoid's equatorial and polar radii. 
%[text] lon: geodetic longitude, rad, from $-\\pi$ to $\\pi$. 
%[text] lat: geodetic latitude, rad, from -PI/2 to PI/2. 
%[text] h: geodetic height, in same units as the variable, a.
%[text] ## note
%[text] NA
%[text] ## references 
%[text] The Astronomical Almanac for the Year 1990, pp. K11-K13. 
%[text] Satellite orbits, pp.188 revisions;
%[text] ## revisions
%[text] 20210427  y.yoshimura
%[text] See also teme2Mod, mod2J2000, obliquity, nutation.
function [lon, lat, h] = geocentric2Geodetic(x, y, z, a, f)
TOLERANCE = 1e-8;

% longitude
lon = atan2(y, x);

r = sqrt(x * x + y * y);
e2 = 2.0 * f - f * f; % eccentricity

lat = atan(z / r);
while true
    p1 = lat;
    s = sin(p1);    
    c = 1.0 / sqrt(1.0 - e2 * s * s);
    lat = atan((z + a * c * e2 * s) / r);
    if (abs(p1 - lat) <= TOLERANCE)
        break;
    end
end

h = r / cos(lat) - a * c;
end

%[appendix]{"version":"1.0"}
%---
