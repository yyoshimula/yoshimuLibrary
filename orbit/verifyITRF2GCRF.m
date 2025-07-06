%[text] # verifying ITRF to GCRF transformation (CIO approcach of IAU-2006/2000 reduction)
%[text] ## cf. example 3-14
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. p.220

% clc
clear

format long

%[text] ## given
rITRF = [-1033.4793830
    7901.2952754
    6380.3565958];

vITRF = [-3.225636520
    -2.872451450
    5.531924446];

yy = 2004;
mm = 4;
dd = 6;
hh = 7;
m = 51;
s = 28.386009;

utc = datetime([yy mm dd hh m s]) %[output:4652f754]
%[text] ## time and EOP
% time
jd = gc2jd(yy, mm, dd, hh, m, s);

leapJD = leapS; % load database
deltaAT = dAT(jd, leapJD);

fname = "EOP_20_C04_one_file_1962-now.txt";
EOP = readEOP(fname);

eopData = eop(yy, mm, dd, EOP.dataAll) %[output:499bb3e7]

%[text] ### time calculation
jdTT = utc2tt(jd, deltaAT) %[output:3325785a]
tTT = jd2jdT(jdTT) % Julian century of TT %[output:2a5d22f4]

%[text] ## yoshimuLibrary function
%[text] `eop`の値がfundameと違うので，結果も少し異なる（はず）
tic
dcm = itrf2gcrf(jd, EOP); %[output:5b3a3056]
toc %[output:1994194b]
rGCRFyoshimura = dcm * rITRF %[output:4e406d9c]

%[text] ## Vallado's
rGCRFvallado = [5102.508959 %[output:group:06bee412] %[output:5cd9b45e]
    6123.011403 %[output:5cd9b45e]
    6378.136925] %[output:group:06bee412] %[output:5cd9b45e]

% Calculate the difference between the two GCRF positions
tic
positionDifference = abs(rGCRFyoshimura - rGCRFvallado) %[output:3254cab1]
toc %[output:39447fcd]
relErr = positionDifference ./ abs(rGCRFvallado) .* 100 %[output:529da8d5]
%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:4652f754]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"utc","value":"   2004\/04\/06 07:51:28\n"}}
%---
%[output:499bb3e7]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"eopData","value":"     mjd: 53101\n      xp: -6.828940068004591e-07\n      yp: 1.617042703835935e-06\n    dUT1: -0.439955200000000\n      dX: -4.702692706762499e-10\n      dY: -4.751174074873453e-10\n     lod: 0.001509500000000\n"}}
%---
%[output:3325785a]
%   data: {"dataType":"textualVariable","outputData":{"name":"jdTT","value":"     2.453101828154746e+06\n"}}
%---
%[output:2a5d22f4]
%   data: {"dataType":"textualVariable","outputData":{"name":"tTT","value":"   0.042623631888994\n"}}
%---
%[output:5b3a3056]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"UT1","value":"   2004\/04\/06 07:51:27\n"}}
%---
%[output:1994194b]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.026764 秒です。\n","truncated":false}}
%---
%[output:4e406d9c]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"3","name":"rGCRFyoshimura","rows":3,"type":"double","value":[["5.102508965158951"],["6.123011400456219"],["6.378136922649801"]]}}
%---
%[output:5cd9b45e]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"3","name":"rGCRFvallado","rows":3,"type":"double","value":[["5.102508959000000"],["6.123011403000000"],["6.378136925000000"]]}}
%---
%[output:3254cab1]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"-5","name":"positionDifference","rows":3,"type":"double","value":[["0.615895169175928"],["0.254378119279863"],["0.235019888350507"]]}}
%---
%[output:39447fcd]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 0.006880 秒です。\n","truncated":false}}
%---
%[output:529da8d5]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"-6","name":"relErr","rows":3,"type":"double","value":[["0.120704377811937"],["0.041544609757746"],["0.036847733297997"]]}}
%---
