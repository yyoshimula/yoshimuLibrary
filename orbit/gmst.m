%[text] # Juliand day to Greenwich mean sidereal time
%[text] `jd`: Julian day, day, nx1 vector
%[text] `GMST`: Greenwich mean sidereal time, rad, nx1 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Jean Meeus, "Astronomical Algorithms, 2nd edition", p.87. Eq.(12.3) or (12.4)
%[text] ## revisions
%[text] 20160419  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gast.
function GMST = gmst(jd)
% arguments
%     jd (:,1) {mustBeNumeric}    
% end
T = (jd - 2451545.0) / 36525.0;

%% Eq.(12.4)
GMST = 280.46061837 + 360.98564736629 .* (jd - 2451545.0) ...
    + 0.000387933 .* T.^2 - T.^3 ./ 38710000.0;

GMST = deg2rad(mod(GMST, 360));

%% Theta0, deg Eq.(12.3) in Astrolib
% SIDEREAL_PER_SOLAR = 1.002737909350795;
% jd0 = floor(jd - 0.5) + 0.5;
% T = (jd0 - 2451545.0) / 36525.0;
% GMST2 = 100.46061837 + T .* 36000.770053608 ...
%     + T.^2 .* 0.000387933 - T.^3 ./ 38710000.0;
% GMST2 = mod(GMST2, 360) / 360;
% GMST2 = 2 * pi * (GMST2 + (jd - jd0)*SIDEREAL_PER_SOLAR);

end

%[appendix]{"version":"1.0"}
%---
