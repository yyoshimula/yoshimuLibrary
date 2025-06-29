%[text] # the mean obliquity of the ecliptic for a given epoch
%[text] `jd`: Julian day, day, nx1 vector
%[text] const: orbital constants
%[text] ## note
%[text] NA
%[text] ## references 
%[text] The Astronomical Almanac for the Year 1984.Astronomical Algorithms, 2nd Editions. p.147, Eq.(22.2)
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also nutation, precession, hms2deg.
function  e = obliquity(jd)
% arguments
%     jd (:,1) {mustBeNumeric}
% end


t = (jd - 2451545.0) ./ 36525.0;

% deg
e = 23.0 + hms2deg(0, 26, 21.448) - t .* hms2deg(0,0,46.8150)...
    -t.^2 .* hms2deg(0, 0, 0.00059) + t.^3 .* hms2deg(0, 0, 0.001813);
% rad
e = deg2rad(e);

end

%[appendix]{"version":"1.0"}
%---
