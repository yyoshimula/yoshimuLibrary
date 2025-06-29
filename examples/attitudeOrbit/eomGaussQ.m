%[text] # Equations of motion for attitude and orbit
%[text] 姿勢・軌道計算のための運動方程式
%[text] ## state variables
%[text] ${\\bf x}=\[a,e,i,\\Omega, w, f, {\\bf q}^T, {\\bf \\omega}^T\]^T$
%[text] `a`: semi-major axis, km 
%[text] `e`: eccentricity 
%[text] `inc`: inclination, rad 
%[text] `raan`: right ascension of ascending node, rad 
%[text] `ome`: argument of perigee, rad
%[text] `f`: true anomaly, rad 
%[text] `q`: quaternion where q(4) is the scalar part 
%[text] `w`: angular rate, rad/s
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210209  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, mainATOM.
function dxdt = eomGaussQ(t_, x_, jdIni, sat, const, EGM, earthVSOP, JB, coefs, ELP)
%[text] ## pre-allocation
% pre-allocation
aRSW = zeros(3,1); % translational acceleration at RSW frame
trq = zeros(3,1); % control torque at body-fixed frame
%[text] ## state variables
jd = jdIni + s2day(t_); % Julian day, day

%[text] ### orbit
a = x_(1); % semi-major axis, km
e = x_(2); % eccentricty
inc = x_(3); % inclination, rad
raan = x_(4); % right ascension of ascending node, rad
ome = x_(5); % argument of perigee, rad
f = x_(6); % true anomaly, rad

p = a * (1 - e^2); % semi-latus rectum
h = sqrt(const.GE * p); % orbital angular momentum
%n = sqrt(const.GE / a^3); % mean motion
u = ome + f; % argument of latitude
r = p / (1.0 + e * cos(f));

[rVec, ~] = oe2rv([a, e, inc, raan, ome ,f], 1, const.GE); % km, km/s
rVec = rVec(:);

%[text] ### attitude
qNorm = norm(x_(7:10));
q1 = x_(7) / qNorm; % quaternion, vector part
q2 = x_(8) / qNorm;
q3 = x_(9) / qNorm;
q4 = x_(10) / qNorm; % scalar part
qv = [q1
    q2
    q3]; % vector part
q = [qv
    q4];
wx = x_(11); % angular rate, rad/s
wy = x_(12);
wz = x_(13);
wVec = [wx
    wy
    wz];

%[text] ### longitude, latitude
rECEF = zyx2dcm(gmst(jd), 0, 0) * rVec; %km, position at ECEF frame, 3x1
lon = atan2(rECEF(2), rECEF(1)); % rad, longitude
lat = atan(rECEF(3) / norm(rECEF(1:2))); % rad, latitude
[geoLon, geoLat, geoH] = geocentric2Geodetic(rVec(1), rVec(2), rVec(3), const.RE, const.fE);
%[text] ## DCM
ijk2RSW = zxz2dcm(raan, inc, u);
ijk2B = q2dcm(4, q'); % IJK to body-fixed frame
ecef2NED = zyx2dcm(lon, -(pi/2 + lat), 0);
ijk2ECEF = zyx2dcm(gmst(jd), 0, 0);
%[text] ## Earth's gravitational potential
dcmNutation = nutationDCM(jd,const);% nutation
dcmPrecession = precessionDCM(const.J2000, jd, const); % precession 
tmpDCM = dcmNutation * dcmPrecession;
rsw2TOD = tmpDCM * ijk2RSW';
tod2RSW = rsw2TOD';

% 3-2-1 Euler angle(lam, -phi, 0)でTOD to PEFなので
lam = atan2(tod2RSW(1,2), tod2RSW(1,1)); % azimuth, TOD frameにおけるR軸のazimuth
phi = -1.0 * asin(-tod2RSW(1,3)); % elevation, TOD frameにおけるR軸のelevation

% sez2TOD = ijk2SEZ(lam, phi)'; % SEZ to IJK(TOD) frame
% lam = lam - gast(jd, const); % (local) longitude at ECEF, rad
% lam = lam + (lam < 0.0) .* 2 * pi;

tmp = egm2008(rVec, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % at Cartesian coordinate
% fSEZ = [-tmp(2); tmp(3); 0]; % South方向とspherical coordinateのphi方向は逆なので-1かける
% fEarth = sez2TOD * fSEZ;
fEarth = tod2RSW * tmp';
aRSW = aRSW + fEarth;
%[text] ## Sun's gravitational force
[aSunIJK, sunIJK] = sunG(jd, rVec', const, earthVSOP);
aSunIJK = aSunIJK(:);
sunIJK = sunIJK(:); % sun position at inertial frame, km, 3x1
aSun = ijk2RSW * aSunIJK; % acceleration by sun
aRSW = aRSW + aSun;
%[text] ## Moon's gravitational force
[aMoon, ~] = moonG(jd, rVec', const, ELP);
aMoon = aMoon(:);
aMoon = ijk2RSW * aMoon; % acceleration by moon
aRSW = aRSW + aMoon;

%[text] ## SRP
sunRelI = sunIJK - rVec; % 3x1 vector, from sat to sun vector at IJK, km
sunDist = norm(sunRelI); 
sunRelB = ijk2B * sunRelI; % km, at body-fixed frame 
sat = srpSimple(sat, sunRelB./norm(sunRelB), sunDist.*10^3, const); % SRP 
% sat = selfShadow(sat, sunRelB./norm(sunRelB)); % calc self-shadowing 
fSRP = sum(sat.force,1)'; % N, 3x1 vector at body-fixed frame 
tSRP = sum(sat.torque,1)'; % Nm, 3x1 vector at body-fixed frame

aRSW = aRSW + ijk2RSW * ijk2B' * fSRP ./ sat.m ./ 10^3; 
trq = trq + tSRP;
%[text] ## air drag
% Cd = 2.2; 
% [~, rho] = jacciaBowman(jd, geoLon, geoLat, geoH, const, JB);
% vRel = [h / e * sin(f) / p 
%     h / r - const.WE * r * ijk2RSW(3,3) 
%     const.WE * r * ijk2RSW(2,3)]; % relative velocity at RSW 
% vRel = ijk2B * ijk2RSW' * vRel; % relative velocity at body-fixed 
% tmp = sat.normal * (vRel ./ norm(vRel)); 
% ampF = -0.5 .* Cd .* sat.area .* tmp ./ sat.m .* rho .* norm(vRel).^2; % scalar, at body-fixed 
% ampF = ampF .* (tmp > 0); % Nx1 
% fAirB = ampF .* (vRel ./ norm(vRel))'; % Nx3 matrix for each facet 
% tAir = cross(sat.pos, fAirB); % Nx3 matrix for each facet 
% tAir = sum(tAir, 1)'; % 3x1 vector at body-fixed
% fAirB = sum(fAirB, 1); % 1x3 vector at body-fixed 
% fAir = ijk2RSW * ijk2B' * fAirB'; % 3x1 air drag at RSW
% 
% aRSW = aRSW + fAir; 
% trq = trq + tAir;
%[text] ## gravitational gradient torque
%[text] earth directional (unit) vector, 3x1 vector at body-fixed frame
nGG = ijk2B * (-rVec./norm(rVec)); 
tGG = 3 * const.GEm / (r*10^3)^3 * cross(nGG, sat.J*nGG); % Nm

trq = trq + tGG;
%[text] ## remanent magnetic torques
bNED = geodeticIGRF(jd, geoLat, geoLon, geoH, coefs)'; % 3x1, nT, at NED frame
bIJK = ijk2ECEF' * ecef2NED' * bNED; % nT< at ECI frame
bB = ijk2B * bIJK .* 10^(-9); % 3x1 vector, T, at body-fixed frame
trq_mag = cross(sat.mag, bB); % Nm

trq = trq + trq_mag;
%[text] ## Equations of motion
%[text] Gauss
Gauss = [2 * a * a / h * (e*sin(f)*aRSW(1) + p / r * aRSW(2))
    p / h * (sin(f) * aRSW(1) + (cos(f) + (e + cos(f)) * r / p) * aRSW(2))
    r * cos(u) / h * aRSW(3)
    r * sin(u) / h / sin(inc) * aRSW(3)
    -p/e/h * (cos(f)*aRSW(1) - sin(f)*(1 + r / p)*aRSW(2)) - cos(inc) * r / h  * sin(u)/sin(inc)*aRSW(3)
    h / r^2 + p / e / h * (cos(f) * aRSW(1) - sin(f) * (1 + r / p) * aRSW(2))];

% rotational motion
rotDynamics = -sat.J^(-1) * cross(wVec, (sat.J * wVec)) + sat.J^(-1) * trq;
qKinematics = qKine(4, q', wVec')';

% do not rewrite the script below
%------->state variables [x, y, z, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz]
dxdt = [Gauss
    qKinematics
    rotDynamics];
   
end

%[appendix]{"version":"1.0"}
%---
