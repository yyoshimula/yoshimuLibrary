%[text] # verifying calculation of moon position 
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
year = 2010;
month = 12;
day = 19;
hour = 0;
min = 0;
s = 0;

mjd = mjuliandate(year, month, day, hour, min, s); % % modified Julian day, UTC
dUT1 = deltaUT1(mjd); % ΔUT1, s
jdUTC = gc2jd(year, month, day, hour, min, s); % JD w.r.t. UTC

jdUT1 = jdUTC + s2day(dUT1); % JD w.r.t. UT1 (= UTC + dUT1)
jdTDT = jdUT1 + s2day(ut2tt(jdUT1)); % JD w.r.t. Terrestrial time (TT = TDT)


%[text] ## moon function by y.yoshimura
ELP = elpConst;

tic
moonY = moon(jdTDT, const, ELP); % moon longitude, latitude, and distance
toc %[output:7a523d06]

%[text] ## SPICE Toolkit
%[text] 光行差補正を行わないとMATLABとほぼ一致
loadSpiceK

date = strcat(num2str(year),'-', num2str(month),'-',num2str(day), 'T',...
    num2str(hour), ':', num2str(min),':', num2str(s));

% ephemeris time (ET = TDB)
et = cspice_str2et( date ); % s

jdSPICE = cspice_unitim(et, 'ET', 'JDTDT') %[output:74cd6fd2]

tic
% LT+Sは光行差補正, noneにもできる
[moonSPICE, ~] = cspice_spkpos('MOON', et, 'J2000', 'LT+S', 'EARTH');
moonSPICE = moonSPICE'; % 1x3
toc %[output:62fafcb0]
%[text] ## Aerospace toolbox
%[text] `planetEphemeris`を使用
tic
moonMATLAB = planetEphemeris(jdTDT, 'Earth', 'Moon')
toc %[output:215053c0]

%[text] ## error, km
err1 = moonSPICE - moonMATLAB
err2 = moonSPICE - moonY
err3 = moonMATLAB - moonY

% error norm
norm(err1) %[output:54e684bf]
norm(err2) %[output:6a777c1d]
norm(err3) %[output:1c425750]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":44}
%---
%[output:7a523d06]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.027534 秒です。\n","truncated":false}}
%---
%[output:74cd6fd2]
%   data: {"dataType":"textualVariable","outputData":{"name":"jdSPICE","value":"     2.455549500766018e+06\n"}}
%---
%[output:62fafcb0]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.002498 秒です。\n","truncated":false}}
%---
%[output:215053c0]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.157059 秒です。\n","truncated":false}}
%---
%[output:54e684bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"   0.011179236218788\n"}}
%---
%[output:6a777c1d]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"   1.796076271559132\n"}}
%---
%[output:1c425750]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"   1.791054971785680\n"}}
%---
