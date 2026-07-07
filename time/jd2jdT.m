%[text] # julian century
%[text] Julian date to Julian centuries elapsed since J2000.0 (JD 2451545.0)
%[text] ## inputs
%[text] `jd`: Julian day, day
%[text] ## outputs
%[text] `T`: Julian centuries since J2000.0
%[text] ## note
%[text] ## references 
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function T = jd2jdT(jd)
% arguments
%     jd (:,1) {mustBeNumeric}
% end

T = (jd - 2451545.0) / 36525;

end



%[appendix]{"version":"1.0"}
%---
