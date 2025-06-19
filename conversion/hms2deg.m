%[text] # hour, min, and sec angles to deg 
%[text] ## input
%[text] hour:
%[text] min: 
%[text] sec. 
%[text] ## output
%[text] out: angle, deg
%[text] 20210419 y.yoshimura 
function out = hms2deg(hour, min, sec)

out = hour .* 15.0 + (min .* 60 + sec) ./ 3600;


end

%[appendix]{"version":"1.0"}
%---
