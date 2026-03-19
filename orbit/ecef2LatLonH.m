%[text] # Converting satellite position to geodetic latitude, longitude, and height
%[text] `r`: satellite position vector@ECEF, km, nx3 matrix
%[text] `lat`: geodetic latitude, rad, nx1 vector
%[text] `lon`: longitude, rad, nx1 vector
%[text] `h`: height above reference ellipsoid
%[text] ## note
%[text] NA
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th ed.," pp.169–172.
%[text] Jean Meeus, "Astronomical Algorithms, 2nd ed.," pp.143-148. 
%[text] ## revisions
%[text] 20210601  y.yoshimura
%[text] See also orbitConst.
function [lat, lon, h] = ecef2LatLonH(r, const)

e2 = 0.00669437999014; % e^2

small = 1e-6; % tolerance
iterMax = 10; % itermation max

n = size(r,1);
rNorm = vecnorm(r,2,2);
rIJ = sqrt(r(:,1).^2 + r(:,2).^2); % nx1 vector
%[text] ## longitude
lon = (r(:,3) ./ abs(r(:,3)) .* pi/2) .* (abs(rIJ) < small) ... % near pole
    + atan2(r(:,2), r(:,1)) .* (abs(rIJ) >= small); % nx1
%[text] ## geodetic latitude
tmp  = asin(r(:,3) ./ rNorm);

for j = 1:n
    oldDelta = tmp(j) + 10.0;  % initialization
    i = 0;
    while((abs(oldDelta - tmp(j)) >= small) && (i < iterMax))
        oldDelta = tmp(j);
        c = const.RE ./ (sqrt(1.0 - e2 * sin(tmp(j))^2));
        
        tmp(j) = atan((r(j,3) + c * e2 * sin(tmp(j))) / rIJ(j));
        i = i + 1;
    end    
end

lat = tmp;
%[text] ## height
C = const.RE ./ (sqrt(1.0 - e2 .* sin(lat).^2));
h = (rIJ ./ cos(lat) - C) .* ((pi * 0.5 - abs(lat)) > deg2rad(1))...
    + (r(:,3) ./ sin(lat) -  C .* (1.0 - e2)) .* ((pi * 0.5 - abs(lat)) <= deg2rad(1)); % near pole

end

%[appendix]{"version":"1.0"}
%---
