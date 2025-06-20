%[text] # (wrapper) calculating atomospheric density using Jaccia–Bowman model 2008
%[text] ## inputs
%[text] `jd`: julian day, day
%[text] `alp: right ascension`geocentricc longitude, rad
%[text] `lat:` geocentric latitude, rad
%[text] h`:` geodetic height, km
%[text] `const:` orbital constant
%[text] `JB:` coefficients
%[text] ## outputs
%[text] `temp`: temperature, 
%[text] `rho`: air density
%[text] ## note
%[text] NA
%[text] ## references 
%[text] ## NA
%[text] ## revisions
%[text] 20221009  y.yoshimura, y.yoshimula@gmail.com
%[text] See also jd2mjd, jd2gc.
function [temp, rho] = jacciaBowman(jd, lon, lat, h, const, JB)
% arguments    
%     jd (1,1) {mustBeNumeric}
%     lon (1,1) {mustBeNumeric}
%     lat (1,1) {mustBeNumeric}
%     h (1,1) {mustBeNumeric}
%     const
%     JB
% end

%[text] ## year, day of the year, and modified Juliana date
[year, month, day, hour, minute, sec] = jd2gc(jd); % by y.yoshimura
MJD = jd2mjd(jd); 

%[text] ## read solar indices
% USE 1 DAY LAG FOR F10 AND S10 FOR JB2008
JD = floor(MJD-1+2400000.5);
i = find(JD == JB.SOLdata(3,:),1,'first');
SOL = JB.SOLdata(:,i);
F10 = SOL(4);
F10B = SOL(5);
S10 = SOL(6);
S10B = SOL(7);

% USE 2 DAY LAG FOR M10 FOR JB2008
SOL = JB.SOLdata(:,i-1);
XM10 = SOL(8);
XM10B = SOL(9);

% USE 5 DAY LAG FOR Y10 FOR JB2008
SOL = JB.SOLdata(:,i-4);
Y10 = SOL(10);
Y10B = SOL(11);

doy = finddays(year,month,day,hour,minute,sec);
i = find(year == JB.DTCdata(1,:) & floor(doy) == JB.DTCdata(2,:),1,'first');
DTC = JB.DTCdata(:,i);
ii = floor(hour)+3;
DSTDTC = DTC(ii);
%[text] ### convert point of interest location (rad, km)
%[text] ### convert longitude to RA
[~, ~, UT1_UTC, ~, ~, ~, ~, ~, TAI_UTC] = IERS(JB.EOPdata,MJD,'l', const);
[~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
[DJMJD0, DATE] = iauCal2jd(year, month, day);
TIME = (60*(60*hour+minute)+sec)/86400;
UTC = DATE+TIME;
TT = UTC+TT_UTC/86400;
TUT = TIME+UT1_UTC/86400;
UT1 = DATE+TUT;
GWRAS = iauGmst06(DJMJD0, UT1, DJMJD0, TT, const); % GMST
XLON = lon;
SAT(1) = mod(GWRAS + XLON, 2*pi);
SAT(2) = lat;
SAT(3) = h;

% SET Sun's right ascension and declination (RADIANS)
% Difference between ephemeris time and universal time
% JD = MJD_UTC+2400000.5;
% [year, month, day, hour, minute, sec] = invjday(JD);
% days = finddays(year, month, day, hour, minute, sec);
% ET_UT = ETminUT(year+days/365.25);
% MJD_ET = MJD_UTC+ET_UT/86400;
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_ET);

MJD_TDB = Mjday_TDB(TT);
[~,~,~,~,~,~,~, ...
 ~,~,~,r_Sun,~] = JPL_Eph_DE430(MJD_TDB, JB);
ra_Sun  = atan2(r_Sun(2), r_Sun(1));
dec_Sun = atan2(r_Sun(3), sqrt(r_Sun(1)^2+r_Sun(2)^2));
SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;

%[text] ## compute density KG/M3 rho
[temp, rho] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC);

end


%[appendix]{"version":"1.0"}
%---
