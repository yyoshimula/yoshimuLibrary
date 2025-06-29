%[text] # mean anomaly usin eccentricity and true anomaly
%[text] `e`: eccentricity
%[text] `f`: ture anomaly, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20221110  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
%[text] See also eAnomaly.
function M = meanAnomaly(e, f)
% arguments
%     e (:,1) {mustBeNumeric}
%     f (:,1) {mustBeNumeric}
% end

% eccentric anomaly
E = eAnomaly(e,f);

M = E - e .* sin(E);

M = mod(M, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
