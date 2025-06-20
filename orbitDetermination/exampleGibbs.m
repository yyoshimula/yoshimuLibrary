%[text] # Gibbs method of orbit determination from three position vectors
%[text] 
%[text] ## references 
%[text] Curtis, Howard D. Orbital Mechanics for Engineering Students. Butterworth-Heinemann, 2013, p231-237.

clc
clear

const = orbitConst;

%[text] ## position vectors
r1 = [-294.32, 4265.1, 5986.7]; 
r2 = [-1365.5, 3637.6, 6346.8];
r3 = [-2940.3, 2473.7, 6555.8];

[v1, v2, v3] = gibbs(r1, r2, r3, const.GE)
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
