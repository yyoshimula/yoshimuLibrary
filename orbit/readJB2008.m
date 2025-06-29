%[text] # reading Jaccia–Bowman 2008 coefficients
%[text] Jaccia–Bowman 2008の係数読み込み
%[text] ## inputs
%[text] nothing
%[text] ## outputs
%[text] `PC`
%[text] `EOP:` Earth orientation parameter
%[text] `SOL:` Space weather data
%[text] `DTC:` geomagnetic storm DTC value
%[text] ## note
%[text] ## references 
%[text] ## revisions
%[text] 20210509  y.yoshimura
%[text] See also jacciaBowman.
function [PC, EOP, SOL, DTC] = readJB2008

load('DE430Coeff.mat'); % DE430Coeff;

PC = DE430Coeff;
%[text] ## Earth orientation parameter (EOP)
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
EOP = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

%[text] ## Space weather data
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOL = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);


%[text] ## Geomagnetic storm DTC value
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTC = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);

fclose(fid);

end

%[appendix]{"version":"1.0"}
%---
