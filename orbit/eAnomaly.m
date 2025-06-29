%[text] # eccentric anomaly using eccentricity and true anomaly
%[text] `e`: eccentricity
%[text] `f`: ture anomaly, rad
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20221110  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
%[text] See also eAnomaly.
function E = eAnomaly(e, f)
% arguments
%     e (:,1) {mustBeNumeric}
%     f (:,1) {mustBeNumeric}
% end
f = mod(f, 2*pi);

% eccentric anomaly
E = 2 .* atan(sqrt((1 - e) ./ (1 + e)) .* tan(f ./ 2));

E = mod(E, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
