%[text] # verifying time calculation
%[text] 時間の計算
clc
clear
cls
%[text] ## condition (UTC)
year = 2003;
month = 12;
day = 19;
hour = 0;
min = 0;
s = 0;
date = strcat(num2str(year),'-', num2str(month),'-',num2str(day), 'T',... %[output:group:5916d3d3] %[output:4200c12b]
    num2str(hour), ':', num2str(min),':', num2str(s)) %[output:group:5916d3d3] %[output:4200c12b]
%[text] ## SPICE
loadSpiceK

% ephemeris time (ET = TDB)
et = cspice_str2et( date ); % s

% ET to JD(Terrestrial dynamical time, TDTに対する）
jdSPICE = cspice_unitim(et, 'ET', 'JDTDT')  %[output:753686ad]

%[text] ## yoshimuLibrary
%[text] $UT1 = UTC+\\Delta UT1$
utc = datetime([year month day hour min s]) %[output:18d4507c]

leapJD = leapS;

fname = "EOP_20_C04_one_file_1962-now.txt";
eopDataAll = readEOP(fname);
eopData = eop(year, month, day, eopDataAll) %[output:1fe35afc]

jd = gc2jd(year, month, day, hour, min, s); % JD w.r.t. UTC
deltaAT = dAT(jd, leapJD);
% jdUT1 = jd + s2day(eopData.dUT1); % JD w.r.t. UT1 (= UTC + dUT1)
jd = utc2tt(jd, deltaAT);
%[text] ## Aerospace Toolbox @ MATLAB
% modified Julian day, UTC
mjd = mjuliandate(year, month, day, hour, min, s) % day %[output:239bd760]
%[text] $\\Delta UT1$
dUT1 = deltaUT1(mjd) % s %[output:48dd7405]

mjdMAT = mjd + s2day(dUT1); % modified JD w.r.t. UT1, day
jdMAT = mjd2jd(mjdMAT) % JD w.r.t. UT1 %[output:4fdb3f89]

jdMATtt = juliandate(year, month, day, hour, min, s) + s2day(deltaAT + 32.184);

tdbMAT = tdbjuliandate([year, month, day, hour, min, s])  %[output:355fe75f]
%[text] ## Comparison
%[text] sで表記
%[text] ### SPICE vs yoshimuLibrary
%[text] 微妙にずれる
day2s(jdSPICE - jd) %[output:3c7647e3]
%[text] ### SPICE vs MATLAB
%[text] MATLABの計算はUT1のままなのでずれる．具体的には↓
%[text] $TT = TAI + 32.184 = UTC + \\Delta AT + 32.184$
day2s(jdSPICE - jdMATtt) %[output:977c8ac1]
day2s(jd - jdMATtt) %[output:249f6c75]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:4200c12b]
%   data: {"dataType":"textualVariable","outputData":{"name":"date","value":"'2003-12-19T0:0:0'"}}
%---
%[output:753686ad]
%   data: {"dataType":"textualVariable","outputData":{"name":"jdSPICE","value":"     2.452992500742870e+06\n"}}
%---
%[output:18d4507c]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"utc","value":"   2003\/12\/19\n"}}
%---
%[output:1fe35afc]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"eopData","value":"     mjd: 52992\n      xp: 3.465157304362297e-07\n      yp: 7.636348772524413e-07\n    dUT1: -0.384231300000000\n      dX: -3.054326190990077e-10\n      dY: -7.417649320975900e-10\n     lod: 7.908000000000000e-04\n"}}
%---
%[output:239bd760]
%   data: {"dataType":"textualVariable","outputData":{"name":"mjd","value":"       52992\n"}}
%---
%[output:48dd7405]
%   data: {"dataType":"textualVariable","outputData":{"name":"dUT1","value":"  -0.384227500000000\n"}}
%---
%[output:4fdb3f89]
%   data: {"dataType":"textualVariable","outputData":{"name":"jdMAT","value":"     2.452992499995553e+06\n"}}
%---
%[output:355fe75f]
%   data: {"dataType":"textualVariable","outputData":{"name":"tdbMAT","value":"     2.452992499999995e+06\n"}}
%---
%[output:3c7647e3]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"     0\n"}}
%---
%[output:977c8ac1]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"     0\n"}}
%---
%[output:249f6c75]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"     0\n"}}
%---
