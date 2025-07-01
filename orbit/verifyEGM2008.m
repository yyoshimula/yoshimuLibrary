%[text] # verifying egm2008 function with Aerospace toolbox
%[text] ## note
%[text] ## references 
%[text] ## revisions
%[text] ©2022  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
clc
clear 
cls

format long
%[text] ## condition
rVec = 8000 .* rand(1,3); % km

%[text] ## read orbital constant and EGM2008 coefficients
const = orbitConst; % orbital constants
EGM.GEODEG = 8; % geoid degree
EGM = readEGM2008(EGM, EGM.GEODEG);
%[text] ## egm2008 function by y.yoshimura
tic
% at ECEF frame described with Cartesian coordinate
acc = egm2008(rVec, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % km/s^2
toc %[output:057c09f0]

% The egm2008 function does not include two-body acceleration, so add the two-body acceleration.
acc = acc + -const.GE / norm(rVec)^3 .* rVec; 
%[text] $\\rm \nm/s^2$  to $\\rm \nkm/s^2$
acc = acc .* 10^3 % m/s^2 %[output:2ec1d6dd]
%[text] ## Aerospace toolbox
%[text] deafultでEGM2008を使う
tic
[gx, gy, gz] = gravitysphericalharmonic(rVec .* 10^3, EGM.GEODEG) %[output:4f2e6e8a] %[output:4d2b1a6f] %[output:0467c6bd]
toc %[output:28809712]

absErr = abs(acc - [gx, gy, gz]) %[output:2a6704b3]
relErr = absErr ./ abs(acc) .* 100 %[output:9b9defc0]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:057c09f0]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.010936 秒です。\n","truncated":false}}
%---
%[output:2ec1d6dd]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"acc","rows":1,"type":"double","value":[["-4.104919010905021","-2.842008551067258","-0.439094230255747"]]}}
%---
%[output:4f2e6e8a]
%   data: {"dataType":"textualVariable","outputData":{"name":"gx","value":"  -4.104919013994521\n"}}
%---
%[output:4d2b1a6f]
%   data: {"dataType":"textualVariable","outputData":{"name":"gy","value":"  -2.842008553206249\n"}}
%---
%[output:0467c6bd]
%   data: {"dataType":"textualVariable","outputData":{"name":"gz","value":"  -0.439094230586224\n"}}
%---
%[output:28809712]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.235778 秒です。\n","truncated":false}}
%---
%[output:2a6704b3]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-8","name":"absErr","rows":1,"type":"double","value":[["0.308950021121746","0.213899076229040","0.033047697911570"]]}}
%---
%[output:9b9defc0]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-7","name":"relErr","rows":1,"type":"double","value":[["0.752633658059993","0.752633471664661","0.752633390156863"]]}}
%---
