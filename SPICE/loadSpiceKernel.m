%[text] # load SPICE kernels
% Clear any existing kernels and prepare for loading new ones
cspice_kclear;

% Define the original directory path for your kernel library
originalDir = '/Users/yyoshimula/Dropbox/MATLAB/libraries/SPICE/mice/kernel/';

% Load a leapseconds file
cspice_furnsh(strcat(originalDir,'earth_assoc_itrf93.tf'));
cspice_furnsh(strcat(originalDir, 'naif0012.tls') );

% Load Earth binary PCK
cspice_furnsh(strcat(originalDir,'earth_200101_990628_predict.bpc'));
cspice_furnsh(strcat(originalDir,'pck00010.tpc'));

% Load planetary ephemeris
cspice_furnsh(strcat(originalDir,'de421.bsp'));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
