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
[EGM.Cnm, EGM.Snm] = readEGM2008('EGM2008_to2190_TideFree.txt', EGM.GEODEG);
%[text] ## egm2008 function by y.yoshimura
tic
% at ECEF frame described with Cartesian coordinate
acc = egm2008(rVec, EGM.GEODEG, EGM.Cnm, EGM.Snm, const); % km/s^2
toc %[output:1685a701]

% The egm2008 function does not include two-body acceleration, so add the two-body acceleration.
acc = acc + -const.GE / norm(rVec)^3 .* rVec; 
%[text] $\\rm \nm/s^2$  to $\\rm \nkm/s^2$
acc = acc .* 10^3 % m/s^2 %[output:0f840046]
%[text] ## Aerospace toolbox
%[text] deafultでEGM2008を使う
tic
[gx, gy, gz] = gravitysphericalharmonic(rVec .* 10^3, EGM.GEODEG) %[output:05d0eb9f] %[output:1d7155c0] %[output:22883e66]
toc %[output:2c906555]

absErr = abs(acc - [gx, gy, gz]) %[output:9da6ec71]
relErr = absErr ./ abs(acc) .* 100 %[output:738e7657]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:1685a701]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.000691 秒です。\n","truncated":false}}
%---
%[output:0f840046]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"acc","rows":1,"type":"double","value":[["-0.110331417751748","-2.623303332966638","-2.889199863131134"]]}}
%---
%[output:05d0eb9f]
%   data: {"dataType":"textualVariable","outputData":{"name":"gx","value":"  -0.110331417834787\n"}}
%---
%[output:1d7155c0]
%   data: {"dataType":"textualVariable","outputData":{"name":"gy","value":"  -2.623303334941023\n"}}
%---
%[output:22883e66]
%   data: {"dataType":"textualVariable","outputData":{"name":"gz","value":"  -2.889199865305642\n"}}
%---
%[output:2c906555]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.007210 秒です。\n","truncated":false}}
%---
%[output:9da6ec71]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-8","name":"absErr","rows":1,"type":"double","value":[["0.008303908949348","0.197438509985659","0.217450812911579"]]}}
%---
%[output:738e7657]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-7","name":"relErr","rows":1,"type":"double","value":[["0.752633213508756","0.752633168663648","0.752633335223544"]]}}
%---
