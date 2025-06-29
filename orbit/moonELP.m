%[text] # Moon's geocentric longitude, latitude, and distance w.r.t. referred to the mean ecliptic and equinox of date.
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] ELPcoeffA, ELPcoeffB: ELP coefficients
%[text] ## outputs
%[text] lon: Moon's geocentric longitude, rad
%[text] lat: Moon's geocentric latitude, rad
%[text] r: Moon's geocentric distance, km
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Jean Meeus, "Astronomical Algorithms, 2nd edition", p. 337.
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sun.
function [lon, lat, r] = moonELP(jd, ELP)
% arguments
%     jd (:,1) {mustBeNumeric}
%     ELP
% end
%[text] ## set coefficients
% for longitude and distance
coeffD = ELP.a(:,1);
coeffM = ELP.a(:,2);
coeffMp = ELP.a(:,3);
coeffF = ELP.a(:,4);
sumL = ELP.a(:,5);
sumR = ELP.a(:,6);

% for latitude
coeffDlat = ELP.b(:,1);
coeffMlat = ELP.b(:,2);
coeffMplat = ELP.b(:,3);
coeffFlat = ELP.b(:,4);
sumB = ELP.b(:,5);
%[text] ## calc
t = ( jd - 2451545.0 ) / 36525.0;

% mean longitude, deg
Lprime = 218.3164477 + 481267.88123421 .* t - 0.0015786 .* t.^2 + t.^3 ./ 538841.0 - t.^4 / 65194000.0;
Lprime = mod(Lprime, 360);

% mean elongation, deg
D = 297.8501921 + t .* 445267.1114034 - t.^2 .* 0.0018819 + t.^3 ./ 545868.0 - t.^4 ./ 113065000.0;
D = mod(D, 360);

% sun's mean anomaly, deg
M = 357.5291092 + t .* 35999.0502909 - t.^2 .* 0.0001536 + t.^3 ./ 24490000.0;
M = mod(M, 360);

% moon's mean anomaly, deg
Mprime = 134.9633964 + t .* 477198.8675055 + t.^2 .* 0.0087414 + t.^3 ./ 69699.0 - t.^4 ./ 14712000.0;
Mprime = mod(Mprime, 360);

% moon's argument of latitude, deg
F = 93.2720950 + t .* 483202.0175233 - t.^2 .* 0.0036539 - t.^3 ./ 3526000.0 + t.^4 ./ 863310000.0;
F = mod(F, 360);

Lprime = deg2rad(Lprime);
D = deg2rad(D);
M = deg2rad(M);
Mprime = deg2rad(Mprime);
F = deg2rad(F);

% additional terms
A = [mod(119.75 + 131.849 * t, 360.0 ), mod(53.09 + 479264.290 * t, 360.0), mod( 313.45 + 481266.484 * t, 360.0 )];
A = deg2rad(A);

% Eq. (47.6)
E = 1.0 - t .* 0.002516 - t.^2 * 0.0000074;

%[text] ## longitude and distance
tmpLon = sin(coeffD .* D + coeffM .* M + coeffMp .* Mprime + coeffF .* F);
tmpR = cos(coeffD .* D + coeffM .* M + coeffMp .* Mprime + coeffF .* F);
%[text] #### $M$or $-M\n$を含むものには$E$を$2M$or$-2M$を含むものは$E^2$をかける
tmp = (abs(coeffM) == 0) .* tmpLon + (abs(coeffM) == 1) .* E .* tmpLon + (abs(coeffM) == 2) .* E.^2 .* tmpLon;
lon = sum(sumL .* tmp, 1);

tmp = (abs(coeffM) == 0) .* tmpR + (abs(coeffM) == 1) .* E .* tmpR + (abs(coeffM) == 2) .* E.^2 .* tmpR;
r = sum(sumR .* tmp, 1);
%[text] ## latitude
tmpLat = sin(coeffDlat .* D + coeffMlat .* M + coeffMplat .* Mprime + coeffFlat .* F);
tmp = (abs(coeffMlat) == 0) .* tmpLat + (abs(coeffMlat) == 1) .* E .* tmpLat + (abs(coeffMlat) == 2) .* E.^2 .* tmpLat;
lat = sum(sumB .* tmp, 1);
%[text] ## additive terms
% deg
lon = lon + 3958 * sin(A(1)) +  1962 * sin(Lprime - F) + 318 * sin(A(2));

lat = lat -2235 * sin(Lprime) + 382 * sin(A(3)) + 175 * sin(A(1) - F) ...
    + 175 * sin(A(1) + F) + 127 * sin(Lprime - Mprime) - 115 * sin(Lprime + Mprime);

lon = Lprime + deg2rad(lon / 10^6); %rad
lat = deg2rad(lat / 10^6); % rad
r = 385000.56 + (r / 10^3); % km

end

%[appendix]{"version":"1.0"}
%---
