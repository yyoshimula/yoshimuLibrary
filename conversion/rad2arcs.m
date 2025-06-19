%[text] # converting arcsecond to radian
%[text] `angle`: angle to be converted, rad
%[text] `rad`: angle converted, arcs
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also arcs2rad.
function arcs = rad2arcs(angle)

arcs = angle .* 3600.0 .* 180.0 ./ pi;

end

%[appendix]{"version":"1.0"}
%---
