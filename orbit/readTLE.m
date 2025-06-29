%[text] # Reading TLE
%[text] ## inputs
%[text] TLE(Two line elements)を読み込み
%[text] `tleName`: Two line elements file's name
%[text] `cons`t: constant parameters for orbit propagation
%[text] spiceFlag: on: if SPICE toolkit is used
%[text]                  off (default): not used
%[text] ## outputs
%[text] `tle:` includes the following values
%[text] `tle.satName`: satellite name
%[text] `tle.satI`D: satellite ID
%[text] `tle.launchYear:` launch year
%[text] `tle.launchNum:` launch number
%[text] `tle.launchPiece:` launchPiece;
%[text] `tle.oe`: mean Orbital elements as
%[text]  `a`: semi-major axis, km
%[text]  `e`: eccentricity
%[text]  `i`: inclination, rad
%[text]  `W`: longitude of the ascending node, rad
%[text]  `w`: argument of perigee, rad
%[text]  `M`: mean anomaly at EPOCH, rad
%[text] `tle.jd:` Julian day of Terrestrial Time (TT) \[day\]
%[text] `tle.n`: mean motion \[rev/day\]Revolution Number at Epoch, s
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20190201  y.yoshimura
%[text] See also orbitConst
function tle = readTLE(tleName, const, tool)
arguments
    tleName {mustBeText}
    const
    tool char {mustBeMember(tool,{'yoshimuLibrary', 'SPICE', 'MATLAB'})} = 'yoshimuLibrary'
end

%[text] ### scan TLE file
data = importdata(tleName, '');
n = size(data,1); % the number of data

tmp = cell2mat(data(1));
if strcmp(tmp(1), '1') % if satellite name is not included
    line1 = cell2mat(data(1:2:n)); % 2行おきに抽出
    line2 = cell2mat(data(2:2:n));
else % if satellite name is included
    line0 = cell2mat(data(1:3:n));
    line1 = cell2mat(data(2:3:n)); % satellite nameがあるときは3行おきに抽出
    line2 = cell2mat(data(3:3:n));

    tle.satName = line0;
end
%[text] ### Information from LINE1
satID = line1(:,3:7);
launchYear = str2num(line1(:,10:11));
launchNum = str2num(line1(:,12:14));
launchPiece = str2num(line1(:,15:17));
year = str2num(line1(:,19:20));
epochDayOfYear = str2num(line1(:,21:32));
tmp = (year > 57) .* 1900 ...
    + (year <= 57) .* 2000;
epochYear = tmp + year; % 4-digit year of double format

if mod(n, 2) == 1 % if object name is included
    tmp = ones((n - 1) / 2, 1);
else
    tmp = ones(n/2, 1);
end

%[text] ### time
% MATLABのdatetime配列
epochDayOfYear = duration(epochDayOfYear*24, 0, 0); % h:m:s
tmpDate = datetime(epochYear, tmp, 0);
ep = tmpDate + epochDayOfYear; % full epoch
jd = gc2jd(ep.Year, ep.Month, ep.Day, ep.Hour, ep.Minute, ep.Second); % julian day @UTC, day

if strcmp(tool, 'yoshimuLibrary') % when yoshimuLibrary is used
    leapJD = leapS; % load database
    deltaAT = dAT(jd, leapJD);
    jd = utc2tt(jd, deltaAT); % JD of terrestrial time

elseif strcmp(tool, 'SPICE') % when SPICE toolkit is used
    loadSpiceK
    % ephemeris time (ET = TDB)
    et = cspice_str2et( date ); % s

    % ET to JD(Terrestrial dynamical time (TDT)に対する）
    jd = cspice_unitim(et, 'ET', 'JDTDT');
    cspice_kclear
    
elseif strcmp(tool, 'MATLAB') && exist('deltaUT1', 'file') == 2 % when Aerospace toolbox is used and if it is installed
    mjd = jd2mjd(jd);
    dUT1 = deltaUT1(mjd); % s
    jd = jd + s2day(dUT1); % JD w.r.t. UT1 (= UTC + dUT1)
    jd = jd + s2day(ut2tt(jd)); % JD w.r.t. Terrestrial time (TT = TDT)

else
    % UT1 to TT
    jd = jd + s2day(ut2tt(jd));
end

%[text] ### Information from LINE2
rpd = str2num(line2(:,53:63));           % rounds per day

oe = [(const.GEday .* (1 ./ rpd ./ 2 ./pi).^2).^(1/3) ...  % semimajor axis
    str2num(strcat(tmp.*'0.', line2(:,27:33)))...      % eccentricity
    str2num(line2(:,9:16)) ...            % [deg] inclination
    str2num(line2(:,18:25))...         % [deg] longitude (or right ascension) of the ascending node
    str2num(line2(:,35:42)) ...            % [deg] argument of periapsis (or perigee)
    str2num(line2(:,44:51))];          % [deg] mean anomaly at epoch

%[text] $i, \\Omega, \\omega, M \\in \[0, 2\\pi\]$
oe(:,3:6) = mod(oe(:,3:6), 360);
oe(:,3:6) = deg2rad(oe(:,3:6)); % nx6 matrix

%[text] ### structure
tle.satID = satID;
tle.launchYear = launchYear;
tle.launchNum = launchNum;
tle.launchPiece = launchPiece;
tle.oe = oe;
tle.jd = jd;
tle.n = rpd; % mean motion [rev/day]

end

%[appendix]{"version":"1.0"}
%---
