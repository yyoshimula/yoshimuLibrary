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
ELP = readELP();

tic
moonY = moon(jdTDT, const, ELP); % moon longitude, latitude, and distance
toc %[output:2c247631]

%[text] ## SPICE Toolkit
%[text] 光行差補正を行わないとMATLABとほぼ一致
loadSpiceKernel

date = strcat(num2str(year),'-', num2str(month),'-',num2str(day), 'T',...
    num2str(hour), ':', num2str(min),':', num2str(s));

% ephemeris time (ET = TDB)
et = cspice_str2et( date ); % s

jdSPICE = cspice_unitim(et, 'ET', 'JDTDT');

tic
% LT+Sは光行差補正, noneにもできる
[moonSPICE, ~] = cspice_spkpos('MOON', et, 'J2000', 'LT+S', 'EARTH');
moonSPICE = moonSPICE'; % 1x3
toc %[output:187f4cd6]
%[text] ## Aerospace toolbox
%[text] `planetEphemeris`を使用
tic
moonMATLAB = planetEphemeris(jdTDT, 'Earth', 'Moon') %[output:7dc86cf3]
toc %[output:2a174648]

%[text] ## error, km
err1 = moonSPICE - moonMATLAB %[output:673fa482]
err2 = moonSPICE - moonY %[output:3bff9328]
err3 = moonMATLAB - moonY %[output:94db0be1]

% error norm
norm(err1) %[output:38a6ca6c]
norm(err2) %[output:328104e8]
norm(err3) %[output:31eb9818]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":44}
%---
%[output:2c247631]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.001319 秒です。\n","truncated":false}}
%---
%[output:187f4cd6]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.000868 秒です。\n","truncated":false}}
%---
%[output:7dc86cf3]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"5","name":"moonMATLAB","rows":1,"type":"double","value":[["2.055243472785064","2.922502596544788","1.488723767646003"]]}}
%---
%[output:2a174648]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.002945 秒です。\n","truncated":false}}
%---
%[output:673fa482]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"err1","rows":1,"type":"double","value":[["11.249488812900381","13.755978224449791","7.172115581546677"]]}}
%---
%[output:3bff9328]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"err2","rows":1,"type":"double","value":[["10.596430880861590","15.400969645357691","6.897531396825798"]]}}
%---
%[output:94db0be1]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"err3","rows":1,"type":"double","value":[["-0.653057932038791","1.644991420907900","-0.274584184720879"]]}}
%---
%[output:38a6ca6c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"  19.162911505775220\n"}}
%---
%[output:328104e8]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"  19.926117353872623\n"}}
%---
%[output:31eb9818]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"   1.791054971785680\n"}}
%---
