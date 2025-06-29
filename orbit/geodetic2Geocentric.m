%[text] # Geocentric position from geodetic latitude, longitude, and height
%[text] `lat`: geodetic latitude, rad, nx1 vector
%[text] `lon`: longitude, rad, nx1 vector
%[text] `h`: height above mean sea level
%[text] RE: earth equatorial radius, km or m
%[text] f: Earth flattening, scalar
%[text] `r`: satellite position vector@ECEF, km or m, nx3 matrix
%[text] ## note
%[text] units of h and RE must be consistent
%[text] 位置ベクトルの単位は，REとhの単位と同じになる．
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210601  y.yoshimura
%[text] See also orbitConst.
function r = geodetic2Geocentric(lat, lon, h, RE, f)
% arguments
%     lat (:,1) {mustBeNumeric}
%     lon (:,1) {mustBeNumeric}
%     h (:,1) {mustBeNumeric}
%     RE (1,1) {mustBeNumeric}
%     f (1,1) {mustBeNumeric}   
% end

tmp1 = RE ./ sqrt(1 - (2*f - f^2) .* sin(lat).^2) + h;

tmp2 = RE .* (1 - f)^2 ./ sqrt(1 - (2*f - f^2) .* sin(lat).^2) + h;

r = [tmp1 .* cos(lat) .* cos(lon), tmp1 .* cos(lat) .* sin(lon), tmp2 .* sin(lat)];


end


%[appendix]{"version":"1.0"}
%---
