%[text] # julian century
%[text] 
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
