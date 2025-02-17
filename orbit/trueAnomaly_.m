% calculating true anomaly and geocentric distance
% inputs
% a: semi-major axis, km, nx1 vector
% e: eccentricity, nx1 vector
% M: mean anomaly, rad, nx1 vector
% outputs
% f: true anomaly, rad, nx1 vector
% r: geocentric distance, km, nx1 vector
% note
% NA
% references
% NA
% revisions
% 20160630  y.yoshimura, y.yoshimula@gmail.com
% See also keplerEq.

function [f, r]= trueAnomaly_(a, e, M)

% eccentric anomaly
E = zeros(length(M), 1);
for i = 1:length(M)
    E(i,1) = keplerEq_(M(i), e(i));
end

r = a .* ( 1 - e .* cos(E));
f = 2 .* atan(sqrt((1 + e) ./ (1 - e)) .* tan(E ./ 2));

f = mod(f, 2*pi);
end
