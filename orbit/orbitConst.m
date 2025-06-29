%[text] # constant parameters for orbit propagation
%[text] 軌道運動に関する定数を構造体`const`に定義
%[text] ## variables
%[text] `const.GE:` gravitational constant of Earth, $\\rm km^3/s^2$
%[text] `const.GEm:` gravitational constant of Earth, $\\rm m^3/s^2$
%[text] `const.GEday:` gravitational constant of Earth, $\\rm km^3/day^2$
%[text] `const.J2`: J2 coefficient
%[text] `const.RE:` Earth radius, km, WGS-84
%[text] `const.REm:` Earth radius, m, WGS-84
%[text] `const.WE:` Earth rotation rate, rad/s, WGS-84
%[text] `const.fE:` earth flattening
%[text] `const.EPS0:` The maen obliquity for J2000.0, rad
%[text] `const.GSday:` gravitaional constant of Sun, $\\rm km^3/day^2$
%[text] `const.GS:` gravitaional constant of Sun, $\\rm km^3/s^2$
%[text] `const.RS:` Sun radius, km
%[text] `const.RSm:` Sun radius, m
%[text] `const.AU:` distance between Sun and Earth, 1 AU, km
%[text] `const.AUm:` distance between Sun and Earth, 1 AU, m
%[text] `const.S0:` Solar constant, $\\rm W/m^2$
%[text] `const.GM:` gravitational constant of moon $\\rm km^3/s^2$
%[text] `const.GMm:` gravitational constant of moon $\\rm m^3/s^2$
%[text] `const.GMday:`  gravitational constant of moon $\\rm km^3/day^2$
%[text] `const.J2000:` Julian day at J2000.0
%[text] `const.c:` light speed, m/s
%[text] ## note
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th edition
%[text] Montenbruck  Oliver  & Eberhard Gill, Satellite Orbits. Springer Science & Business Media  2012.
%[text] ## revisions
%[text] 20210419  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function const = orbitConst()
%[text] ## Earth
const.GE = 398600.4415;  % km^3/s^2
const.GEm = const.GE * 10^3 * 10^3 * 10^3; % m^3/s^2
const.GEday = const.GE *  86400 * 86400;  % GE, km^3/day^2

const.J2 = 0.00108263;
const.RE = 6378.137; % Earth radius, km, WGS-84
const.REm = const.RE * 10^3; % Earth radius, km, WGS-84
const.WE = 7.292115146706979 * 10^(-5); % Earth rotation rate, rad/s, WGS-84

const.fE = 1.0 / 298.257; % earth flattening
%[text] ## Sun
% The maen obliquity for J2000.0, rad
const.EPS0 = arcs2rad(84381.448);

const.GSday = 990693056236769280000.0;  % km^3/day^2
const.GS = const.GSday / 86400 / 86400; % km^3/s^2
const.RS = 696000.0; % Sun radius, km
const.RSm = const.RS * 10^3; % Sun radius, km
const.AU = 149597870.691; % 1 AU, km
const.AUm = const.AU * 10^3; % 1 AU, m

const.S0 = 1357; % Solar constant, W/m^2
%[text] ## moon
const.GM = 4902.801;     % km^3/s^2
const.GMm = const.GM * 10^3 * 10^3 * 10^3;     % m^3/s^2
const.GMday = 36599213352960.0;     % km^3/day^2
%[text] ## other
const.J2000 = 2451545.0; % Julian day
const.c = 299792458; % light speed, m/s

% save('orbit_const.mat')

end

%[appendix]{"version":"1.0"}
%---
