%[text] # Julian date to modified Julian date
%[text] ## inputs
%[text] `jd:` Julian day, day
%[text] ## outputs
%[text] `mjd:` modified Julian day, day
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230202  y.yoshimura
%[text] See also orbitConst.
function mjd = jd2mjd(jd)
% arguments
%     jd (:,1) {mustBeNumeric}
% end

mjd = jd - 2400000.5;

end

%[appendix]{"version":"1.0"}
%---
