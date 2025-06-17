%[text] # calculating the non-spherical Earth's gravitational attraction with EGM2008 geopotential
%[text] EGM2008を用いた高次の地球重力加速度
%[text] ## input
%[text] `rVec`: object's position at ECEF(Earth-centered Earth-fixed) frame, 1x3 vector, km 
%[text] `deg`: Degree and order required (up to degree and order 20). 
%[text] `Cnm`, `Snm`: Earth's gravitational coefficients of EGM2008
%[text] `const`: orbital constants
%[text] ## output
%[text] `a`: perturbing accelerations in ECEF frame (Cartesian coordinates), 1x3 vector, $\\rm km/s^2\n$
%[text] ## note
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th edition, pp.538–550.
%[text] Montenbruck  Oliver  & Eberhard Gill, Satellite Orbits. Springer Science & Business Media  2012., pp.61–67
%[text] ## how to use
%[text] cf. verifyEGM2008.m
%[text] ```matlabCodeExample
%[text] EGM.GEODEG = 8; % choose the geoid degree
%[text] 
%[text] % read orbital constant and EGM2008 coefficients
%[text] const = orbitConst;
%[text] [EGM.Cnm, EGM.Snm] = readEGM2008('EGM2008_to2190_TideFree.txt', EGM.GEODEG);
%[text] acc = egm2008(rVec, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % at ECEF frame and Cartesian coordinate
%[text] ```
%[text] ## revisions
%[text] 20230119 arguments changed, y.yoshimura
%[text] 20210419  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readEGM2008, orbitConst, verifyEGM2008.
function a = egm2008(rVec, deg, Cnm, Snm, const)
arguments
    rVec (1,3) {mustBeNumeric}
    deg (1,1) {mustBeNumeric}
    Cnm
    Snm
    const
end

x = rVec(1);
y = rVec(2);
z = rVec(3);
r = norm(rVec);

phi = asin(rVec(:,3) ./ r); %geocentric latitude
lam = atan2(rVec(:,2), rVec(:,1)); % longitude
sp = sin(phi);

%[text] ## associated Legendre polynomials, P 
%[text] MATLABのlegendre function
%[text] $P\_{n,m}(x)=(-1)^{m}(1-x^{2})^{m/2}\\frac{\\mathrm{d}^{m}}{\\mathrm{d}x^{m}}P\_{n}(x)$
P = zeros(deg+1, deg+1); % matlabはindexが1始まりなので+1

for i = 1:deg
    tmp = 0:i; % malabのlegendre functionは(-1)^mが付くため
    P(i+1,1:i+1) = (-1).^tmp .* legendre(i, sp)';
end
P(1,:) = 0;
P(2,:) = 0; % 2次以上（degree >=2）しか使わないので0にする．

%[text] ## accelerations (except for two-body acceleration)
%[text] ### partial derivatives of potential
% Clm: (n+1) x (m+1) matrix
% Slm: (n+1) x (m+1) matrix
%[text] $\\frac{\\partial U}{{\\partial r}$
% dU/dr
tmp = lam * (0:deg);
tmpC = repmat(cos(tmp), deg+1, 1); % cos(m*lam)が並んだ行列, 列方向にm, (deg+1) x (deg+1) matrix
tmpS = repmat(sin(tmp), deg+1, 1); % cos(m*lam)が並んだ行列

n = 0:deg; % degree of geopotential
n = n(:); % column vector
tmpCoef = (const.RE ./ r).^n .* (n + 1) .* P .* (Cnm .* tmpC + Snm .* tmpS); %row: degree, column: order
dudr = -const.GE / r^2 * sum(tmpCoef, 'all');

%[text] $\\frac{\\partial U}{{\\partial \\phi}$
% dU/dphi
P2 = [P(:,2:end) zeros(deg+1,1)]; % P(n,m+1)として使う, (deg+1) x (deg+1) matrix
m = repmat(0:deg, deg+1, 1); % order of geopotential, (deg+1) x (deg+1) matrix
tmpCoef = (const.RE / r).^n .* (P2 - tan(phi) * m .* P) .* (Cnm .* tmpC + Snm .* tmpS);
dudphi = const.GE / r * sum(tmpCoef, 'all');
%[text] $\\frac{\\partial U}{\\partial \\lambda}$
% dU/dlambda
tmpCoef = (const.RE / r).^n .* m .* P .* (-Cnm .* tmpS + Snm .* tmpC);
dudlam = const.GE / r * sum(tmpCoef, 'all');

%[text] ### accelerations w.r.t. Earth-fixed frame (except for two-body acceleration)
% acceleration
a = [(1 / r * dudr - z / r^2 / sqrt(x^2 + y^2) * dudphi) * x - y / (x^2 + y^2) * dudlam, ...
    (1 / r * dudr - z / r^2 / sqrt(x^2 + y^2) * dudphi) * y + x / (x^2 + y^2) * dudlam, ...
    z / r * dudr + sqrt(x^2 + y^2) / r^2 * dudphi];

end

%[appendix]{"version":"1.0"}
%---
