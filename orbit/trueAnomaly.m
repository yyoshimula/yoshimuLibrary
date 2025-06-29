%[text] # calculating true anomaly and geocentric distance
%[text] ## inputs
%[text] `a:` semi-major axis, km, nx1 vector
%[text] `e:` eccentricity, nx1 vector
%[text] `M`: mean anomaly, rad, nx1 vector
%[text] ## outputs
%[text] `f`: true anomaly, rad, nx1 vector 
%[text] `r`: geocentric distance, km, nx1 vector 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20160630  y.yoshimura, y.yoshimula@gmail.com
%[text] See also keplerEq.
function [f, r]= trueAnomaly(a, e, M)
% arguments
%     a (:,1) {mustBeNumeric}
%     e (:,1) {mustBeNumeric}
%     M (:,1) {mustBeNumeric}
% end

% eccentric anomaly
E = zeros(length(M), 1);

if isscalar(e)
    e = repmat(e, length(M), 1);
end

for i = 1:length(M)
    E(i,1) = keplerEq(M(i), e(i));
end

r = a .* ( 1 - e .* cos(E));
f = 2 .* atan(sqrt((1 + e) ./ (1 - e)) .* tan(E ./ 2));

f = mod(f, 2*pi);
end

%[appendix]{"version":"1.0"}
%---
