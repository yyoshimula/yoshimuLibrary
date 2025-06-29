%[text] # modified Julian date to Julian date
%[text] ## inputs
%[text] `mjd: modified` Julian day, day
%[text] ## outputs
%[text] `jd:` Julian day, day
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230202  y.yoshimura
%[text] See also orbitConst.
function jd = mjd2jd(mjd)
% arguments
%     mjd (:,1) {mustBeNumeric}
% end

jd = mjd + 2400000.5;

end

%[appendix]{"version":"1.0"}
%---
