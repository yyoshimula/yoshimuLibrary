%[text] # Orbit determination by double-r method, not yet
%[text] ## note
%[text] NA
%[text] ## references
%[text] Vallado, D.A., & Wayne D. McClain. Fundamentals of Astrodynamics and Applications. 4th edition, Springer Science & Business Media, 2001. pp447,
%[text] ## example7-2
% clc
clear

format long

const = orbitConst;
%[text] ## kernels for SPICE toolkit
% Load a leapseconds file
% cspice_furnsh('earth_assoc_itrf93.tf');
% cspice_furnsh( 'naif0012.tls' )
% % Load Earth binary PCK
% cspice_furnsh('earth_000101_230228_221205.bpc');
% cspice_furnsh('pck00010.tpc')
% % Load planetary ephemeris
% cspice_furnsh('de421.bsp');

%[text] ## observer and observations
lat = deg2rad(40); % latitude
lon = deg2rad(-110); % longitude
alt = 2000; % m

% right ascension
alp = [0.939913
    45.025748
    67.886655];
alp = deg2rad(alp);

% declination
dlt = [18.667717
    35.664741
    36.996583];
dlt = deg2rad(dlt);

% line of sight direction vectors
[L(:,1), L(:,2), L(:,3)] = sph2cart(alp, dlt, 1);

% observation epoch
dUT1 = -0.609641; % s
dAT = 35; % s
% Earth orientation parameters
xp = arcs2rad(0.137495);
yp = arcs2rad(0.342416);

%[text] ## epoch
% UTC
year_ = 2012 * ones(3,1);
month_ = 8 * ones(3,1);
day_ = 20 * ones(3,1);
hour_ = [11; 11; 11];
min_ = [40; 48; 52];
sec_ = [28; 28; 28];

% Julian day
jd = gc2jd(year_, month_, day_, hour_, min_, sec_)

jdUT1 = jd + s2day(dUT1);
jdTAI = jd + s2day(dAT);
jdTT = jdTAI + s2day(32.184);


UTC1 = [2012 08 20 11 40 28];
%[text] ### FK5-based
% local sidereal time
% theta = gast(jdUT1, const) + deg2rad(lon)
% 
% % observer position at ECEF frame, km
obsECEF = geodetic2Geocentric(lat, lon, alt/1000, const.RE, const.fE);
% 
% qI2TOD = earthNutationPrecessionQ(const.J2000, jdUT1, 4, const);
% qTOD2PEF = [zeros(size(jd,1), 2), sin(theta./2), cos(theta./2)];
% qI2PEF = qMult(4, 1, qTOD2PEF, qI2TOD);
% 
% W = pef2itrf(xp, yp);
% qI2ITRF = qMult(4, 1, repmat(dcm2q(4, W), size(jd,1),1), qI2PEF);
% 
% % FK5でpolar motion等を考慮してないのでずれる？
% obsECI = qRotation(4, repmat(obsECEF, size(jd,1),1), qInv(4, qI2ITRF))

%[text] ### IAU2000-based using yoshimuLibrary
fname = "EOP_20_C04_one_file_1962-now.txt";

leapJD = leapS;
iau06 = readIAU06;
eopDataAll = readEOP(fname);
dcmITRF2GCRF = itrf2gcrf(datetime(UTC1), leapJD, iau06, eopDataAll);

obsECI = dcmITRF2GCRF * obsECEF';
obsECI = obsECI'
%[text] #### use MATLAB function
r= ecef2eci(UTC1, obsECEF)


dcmMATLAB = dcmeci2ecef('IAU-2000/2006', UTC1, dAT, dUT1, [xp, yp]);

dcmMATLAB'*obsECEF'

% values on fundame
obsECItrue = [4054.881, 2748.195, 4074.237 %[output:group:8a843882]
    3956.224, 2888.232, 4074.364
    3905.073, 2956.935, 4074.430] %[output:group:8a843882]

%[text] ## SPICE
% date1 = 'Aug 20 2012 11:40:28';
% date2 = 'Aug 20 2012 11:48:28';
% date3 = 'Aug 20 2012 11:52:28';
% 
% et1 = cspice_str2et( date1 );
% jd1 = cspice_unitim(et1, 'ET', 'JDTDT');
% et2 = cspice_str2et( date2 );
% jd2 = cspice_unitim(et2, 'ET', 'JDTDT');
% et3 = cspice_str2et( date3 );
% jd3 = cspice_unitim(et3, 'ET', 'JDTDT');
% 
% dcm1 = cspice_pxform('J2000', 'ITRF93', et1);
% dcm2 = cspice_pxform('J2000', 'ITRF93', et2);
% dcm3 = cspice_pxform('J2000', 'ITRF93', et3);
% 
% [dcm1' * obsECEF', dcm2' * obsECEF', dcm3' * obsECEF']

%[text] ## Gauss method
% rRange = 0:20000;
% [r2, v2] = gauss(jd, [alp, dlt], obsECI, rRange, const.GE)
%[text] ## double-r method
r2True = [6356.486034, 5290.5322578, 6511.396979];
v2True = [-4.172948, 4.776550, 1.720271];

aTrue = 12246.023;
eTrue = 0.2;
iTrue = 40; % deg
raanTrue = 330; % deg

[r2, v2] = doubleR(jdUT1, [alp, dlt], obsECI, const.GE, const.RE); %[output:8deb1fd3]
[a, e, inc, raan, w, f, M, lon] = rv2oe(r2, v2, const.GE);

r2 
norm(r2)

rad2deg(inc)
rad2deg(raan)
rad2deg(w)
rad2deg(f)
rad2deg(M)


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40.4}
%---
%[output:8deb1fd3]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('dot')\" style=\"font-weight:bold\">dot<\/a>\nA と B は同じサイズでなければなりません。\n\nエラー: <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('doubleR', '\/Users\/yyoshimula\/Dropbox\/MATLAB\/MATLABdrive\/yoshimuLibrary\/orbitDetermination\/doubleR.mlx', 17)\" style=\"font-weight:bold\">doubleR<\/a> (<a href=\"matlab: opentoline('\/Users\/yyoshimula\/Dropbox\/MATLAB\/MATLABdrive\/yoshimuLibrary\/orbitDetermination\/doubleR.mlx',17,0)\">行 17<\/a>)\nc = 2 .* dot(L, rObs, 2);"}}
%---
