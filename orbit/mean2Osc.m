%[text] # Converting mean orbital elements to osculating orbital elements
%[text] 平均軌道要素を接触軌道要素へ変換
%[text] ## inputs
%[text] `n`: mean motion, rev/day
%[text] `e`: eccentricity
%[text] `i`: inclination, rad
%[text] `Ome`: arugment of ascending node, rad
%[text] `w`: argument of perihelion, rad
%[text] `f`: true anomaly, rad
%[text] `M`: mean anomaly, rad
%[text] `const`:  constant parameters for orbit propagation
%[text] ## outputs
%[text] osculating orbital elements: $\[a\_{\\rm osc},e\_{\\rm osc},i\_{\\rm osc}, \\Omega\_{\\rm osc},\\omega\_{\\rm osc},f\_{\\rm osc},M\_{\\rm osc},r\_{\\rm osc},\\dot{r}\_{\\rm osc}, p\_{\\rm osc}, u\_{\\rm osc}\]$
%[text] `osc`: \[aOsc,  eOsc, iOsc, OmeOsc, wOsc, fOsc, MOsc, rOsc, drOsc, pOsc, uOsc\]
%[text] ## note
%[text] NA
%[text] ## references 
%[text] \[1\] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th edition,  pp.708-709.
%[text] ## revisions
%[text] 20210531  y.yoshimura
%[text] See also orbitConst, trueAnomaly.
function osc = mean2Osc(n, e, i, Ome, w, M, const)
%[text] ### mean values
nSec = 2 * pi .* n ./ (24 * 60 * 60); % mean motion, rad/s
a = (const.GE ./ nSec.^2).^(1/3); % semi-major axis, km
q = a .* (1 - e);
p = q .* (1 + e);

[f, ~] = trueAnomaly(a, e, M);
u = w + f; % argument of latitude

r = p ./ (1 + e .* cos(f));
drMean = sqrt(const.GE ./ p) .* e .* sin(f); % radial velocity

coef = const.J2 .* const.RE.^2;

%[text] ### short-period variations
%[text] $\\Delta i\_{SP}, \\Delta p\_{SP}, \\Delta \\Omega\_{SP}$ in \[1\]
%[text] ### inclination, semiparameter, and RAAN
tmp = 3.0 .* cos(2.0 .* u) + 3.0 .* e .* cos(2.0 .* w + f) ...
    + e .* cos(2.0 .* w + 3.0 .* f);

di = tmp .* coef .* sin(i) .* cos(i) / 4 ./ p.^2;
dp = tmp .* coef * sin(i).^2 / 2 ./ p;

dOme = 6.0 * (f - M + e .* sin(f)) - 3.0 .* sin(2.0 * u) ...
    -3.0 * e * sin(2.0 * w + f) - e * sin(2.0 * w + 3.0 * f);
dOme = -dOme .* coef .* cos(i) ./ 4 ./ p.^2;
%[text] ### position, radial velocity, and argument of latitude
%[text] $\\Delta r\_{SP}, \\Delta \\dot{r}\_{SP}, \\Delta u\_{SP}$ in \[1\]
dr = (3 * cos(i).^2 - 1) .* (2*sqrt(1-e.^2)/(1+e.*cos(f)) ...
    + e.*cos(f) ./ (1+sqrt(1-e^2)) + 1) ...
    - sin(i).^2 .* cos(2.*u);
dr = -dr .* coef ./ 4 ./ p;

ddr = (3 * cos(i).^2 - 1) .* e .* sin(f) .* (sqrt(1-e.^2) ...
    + (1+e.*cos(f)).^2 ./ (1+sqrt(1 - e.^2)) ...    
    - 2 .* sin(i).^2 .* (1 - e.*cos(f)).^2 .* sin(2.*u));
ddr = ddr .* coef .* sqrt(const.GE) ./ 4 ./ p^(5/2);

du = (6 - 30.*cos(i).^2).*(f - M) + 4.*e.*sin(f).*((1 - 6*cos(i).^2) ...
    + (1-3.*cos(i).^2) ./ (1 + sqrt(1-e.^2))) ...
    + (1 - 3.*cos(i).^2) ./ (1+sqrt(1-e.^2)) .* e^2 .* sin(2*f) ...
    + (5 * cos(i).^2 - 2) .*2 * e * sin(f+2*w)...
    +(7 * cos(i).^2 - 1) .* sin(2*u) + 2.*cos(i).^2 * e * sin(3*f+2*w);
du = du * coef / 8 / p.^2;

%[text] ### correct mean values
rOsc = r + dr;
drOsc = drMean + ddr;
pOsc = p + dp;

A = pOsc ./ rOsc - 1;
B = sqrt(pOsc /const.GE) .* drOsc;
eOsc = sqrt(A^2 + B^2);
aOsc = pOsc / (1 - eOsc^2);
iOsc = i + di;
OmeOsc = Ome + dOme;
uOsc = u + du;
fOsc = atan2(B, A);
wOsc = uOsc - fOsc;

% eccentric anomaly
ea = 2.0 * atan(sqrt((1.0 - eOsc) / (1.0 + eOsc)) * tan(0.5 * fOsc));

% mean anomaly
MOsc = ea - eOsc * sin(ea);

%[text] $\\Omega, \\omega, f, M \\in \[0, 2\\pi\]$
OmeOsc = mod(OmeOsc, 2*pi);
wOsc = mod(wOsc, 2*pi);
fOsc = mod(fOsc, 2*pi);
MOsc = mod(MOsc, 2*pi);

osc = [aOsc,  eOsc, iOsc, OmeOsc, wOsc, fOsc, MOsc, rOsc, drOsc, pOsc, uOsc];

end

%[appendix]{"version":"1.0"}
%---
