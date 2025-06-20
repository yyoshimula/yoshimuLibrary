%[text] # inverse calculation of apparent magnitude of object with distance d \[m\]
%[text] `m`: relative magnitude of light curves w.r.t. Sun, Nx1 vector
%[text] d: distance between object and observer, m, Nx1 vector
%[text] $f\_{\\rm obs}\n$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] for the definitoin of $f\_{\\rm obs}\n$: `Light Curve Approximation Using an Attitude Model of Solar Sail Spacecraft'       
%[text] ## revisions
%[text] 20200430  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readSC, lcSimple, lcAS, lcCT, orbitConst.
function fobs = magInv(mApp, d)

msun = -26.7; % apparent magnitude of Sun

fobs =  d.^2 .* 10.^((mApp - msun) ./ (-2.5));

end

%[appendix]{"version":"1.0"}
%---
