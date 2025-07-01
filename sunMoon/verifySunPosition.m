%[text] # verifying calculation of Sun position 
%[text] ## note
%[text] ## references 
%[text] ## revisions
%[text] ©2022  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
clc
clear
cls

format long
%[text] ## condition
const = orbitConst;

% UTC
yy = 2010;
mm = 12;
dd = 19;
hh = 0;
m = 0;
s = 0;

utc = datetime([yy mm dd hh m s]) % datetime format %[output:126e6871]

jdUTC = gc2jd(yy, mm, dd, hh, m, s); % Julian day

leapJD = leapS; % load database
deltaAT = dAT(jdUTC, leapJD);

fname = "EOP_20_C04_one_file_1962-now.txt";
eopDataAll = readEOP(fname);
eopData = eop(yy, mm, dd, eopDataAll) %[output:82e860ae]

%[text] ### UT!, TAI and TT
UT1 = utc + seconds(eopData.dUT1) %[output:982e861c]
TAI = utc + seconds(deltaAT) %[output:6b0c4632]
TT = TAI + seconds(32.184) %[output:2e7a53ab]
%[text] ## sun function by y.yoshimura
%[text] MATLABと150 kmくらいの差異
earthVSOP = vsopConst;
jdTT = gc2jd(year(TT), month(TT), day(TT), hour(TT), minute(TT), second(TT));

tic
sunY = sun(jdTT, const, earthVSOP); % longitude, latitude, and distance
toc %[output:87e56b5c]

%[text] ## SPICE Toolkit
%[text] 光行差補正を行わないとMATLABとほぼ一致
loadSpiceKernel

date = strcat(num2str(yy),'-', num2str(mm),'-',num2str(dd), 'T',...
    num2str(hh), ':', num2str(m),':', num2str(s));

% ephemeris time (ET = TDB)
et = cspice_str2et( date ); % s

jdSPICE = cspice_unitim(et, 'ET', 'JDTDT') %[output:53aab568]

tic
% LT+Sは光行差補正, noneにもできる
[sunSPICE, ~] = cspice_spkpos('SUN', et, 'J2000', 'none', 'EARTH');
sunSPICE = sunSPICE'; % 1x3
toc %[output:03455fd4]
%[text] ## MATLAB Aerospace toolbox
%[text] `planetEphemeris`を使用（アドオンでephemeris dataを追加する必要あり）
tic
mjdUTC = mjuliandate(yy, mm, dd, hh, m, s); % % modified Julian day, UTC
dUT1 = deltaUT1(mjdUTC); % ΔUT1, s
jdUT1 = jdUTC + s2day(dUT1); % JD w.r.t. UT1 (= UTC + dUT1)
jdTT = jdUT1 + s2day(ut2tt(jdUT1)); % JD w.r.t. Terrestrial time (TT = TDT)

sunMATLAB = planetEphemeris(jdTT, 'Earth', 'Sun') %[output:012ab026]
toc %[output:0de30bb7]

%[text] ## error, km
err1 = sunSPICE - sunMATLAB %[output:21bb003e]
err2 = sunSPICE - sunY %[output:2ecdb84a]
err3 = sunMATLAB - sunY %[output:7bd7038b]

% error norm
norm(err1) %[output:181da97e]
norm(err2) %[output:79829b5a]
norm(err3) %[output:02a18020]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":44}
%---
%[output:126e6871]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"utc","value":"   2010\/12\/19\n"}}
%---
%[output:82e860ae]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"eopData","value":"     mjd: 55549\n      xp: 7.908668617203638e-07\n      yp: 1.055691486889637e-06\n    dUT1: -0.132660200000000\n      dX: 1.076286372063170e-09\n      dY: -9.162978572970231e-10\n     lod: 1.911000000000000e-04\n"}}
%---
%[output:982e861c]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"UT1","value":"   2010\/12\/18 23:59:59\n"}}
%---
%[output:6b0c4632]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"TAI","value":"   2010\/12\/19 00:00:34\n"}}
%---
%[output:2e7a53ab]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"TT","value":"   2010\/12\/19 00:01:06\n"}}
%---
%[output:87e56b5c]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.003482 秒です。\n","truncated":false}}
%---
%[output:53aab568]
%   data: {"dataType":"textualVariable","outputData":{"name":"jdSPICE","value":"     2.455549500766018e+06\n"}}
%---
%[output:03455fd4]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.000980 秒です。\n","truncated":false}}
%---
%[output:012ab026]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"8","name":"sunMATLAB","rows":1,"type":"double","value":[["-0.081914341369389","-1.348400955645031","-0.584557174683310"]]}}
%---
%[output:0de30bb7]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.039325 秒です。\n","truncated":false}}
%---
%[output:21bb003e]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"err1","rows":1,"type":"double","value":[["-0.655843310058117","0.488258361816406","-1.065824441611767"]]}}
%---
%[output:2ecdb84a]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"2","name":"err2","rows":1,"type":"double","value":[["0.655191783662885","-0.575454847514629","1.260287463739514"]]}}
%---
%[output:7bd7038b]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"2","name":"err3","rows":1,"type":"double","value":[["0.661750216763467","-0.580337431132793","1.270945708155632"]]}}
%---
%[output:181da97e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"   1.343319923014865\n"}}
%---
%[output:79829b5a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"     1.532562901210393e+02\n"}}
%---
%[output:02a18020]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"     1.545965030794528e+02\n"}}
%---
