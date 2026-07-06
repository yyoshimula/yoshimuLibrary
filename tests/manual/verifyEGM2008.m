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

% The egm2008 function does not include two-body acceleration, and add the two-body acceleration.
acc = acc - const.GE / norm(rVec)^3 .* rVec; 
%[text] $\\rm \nkm/s^2$ to $\\rm \nm/s^2$ 
acc = acc .* 10^3 % m/s^2 %[output:2ec1d6dd]
%[text] ## Aerospace toolbox
%[text] deafultでEGM2008を使う
tic
[gx, gy, gz] = gravitysphericalharmonic(rVec .* 10^3, EGM.GEODEG) %[output:40923ec8] %[output:662798e7] %[output:359e2c15] %[output:9a8b0163]
toc %[output:115af6b5]

absErr = abs(acc - [gx, gy, gz]) %[output:019e433a]
relErr = absErr ./ abs(acc) .* 100 %[output:63a34f67]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:057c09f0]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.001810 秒です。\n","truncated":false}}
%---
%[output:2ec1d6dd]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"acc","rows":1,"type":"double","value":[["-1.899698869041642","-24.422484332638849","-21.493889637293730"]]}}
%---
%[output:40923ec8]
%   data: {"dataType":"warning","outputData":{"text":"警告: Radial position is less than equatorial radius of planetary model, 6378137."}}
%---
%[output:662798e7]
%   data: {"dataType":"textualVariable","outputData":{"name":"gx","value":"  -1.899698870471420\n"}}
%---
%[output:359e2c15]
%   data: {"dataType":"textualVariable","outputData":{"name":"gy","value":" -24.422484351020046\n"}}
%---
%[output:9a8b0163]
%   data: {"dataType":"textualVariable","outputData":{"name":"gz","value":" -21.493889653470763\n"}}
%---
%[output:115af6b5]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.003353 秒です。\n","truncated":false}}
%---
%[output:019e433a]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-7","name":"absErr","rows":1,"type":"double","value":[["0.014297778516692","0.183811970089209","0.161770330464606"]]}}
%---
%[output:63a34f67]
%   data: {"dataType":"matrix","outputData":{"columns":3,"exponent":"-7","name":"relErr","rows":1,"type":"double","value":[["0.752633943710507","0.752634202096951","0.752634042485825"]]}}
%---
