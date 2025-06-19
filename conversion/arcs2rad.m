%[text] # converting arcsecond to radian
%[text] `angle`: angle to be converted, arcs
%[text] `rad`: angle converted, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function rad = arcs2rad(angle)

rad = angle .* pi ./ 3600 ./ 180.0;

end

%[appendix]{"version":"1.0"}
%---
