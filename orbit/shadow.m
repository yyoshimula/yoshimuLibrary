%[text] # Earth shadow function
%[text] `satI`: satellite position@inertial frame, arbitrary unit, nx3
%[text] `sunI`: sun position@inertial frame, nx3
%[text] `rS`: Sun radius, scalar
%[text] `rE`: Earth radius, scalar
%[text] nu: 1: sunlit, 0: umbra (eclipse)
%[text] ## note
%[text] all variables (`satI, sunI, rS, rE`) must have the same unit
%[text] 変数の単位は統一すること．（統一していればmでもkmでもok）
%[text] ## references 
%[text] Satellite Orbits, Montenbruck, Gill, p.81
%[text] ## revisions
%[text] 20231128 y.yoshimura: 引数を変更
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also vsopConst.
function  nu = shadow(satI, sunI, rS, rE)
% arguments
%     satI (:,3) {mustBeNumeric}
%     sunI (:,3) {mustBeNumeric}
%     rS (1,1) {mustBeNumeric}
%     rE (1,1) {mustBeNumeric}    
% end

sunRel = sunI - satI;

% apparent Sun disk radius
a = asin(rS ./ vecnorm(sunRel,2,2)); % nx1 vector

% apparent Earth disk radius
b = asin(rE ./ vecnorm(satI, 2, 2));

% apparent separation of the two body centers
c = acos(dot(-satI, sunRel,2) ./ vecnorm(satI,2,2) ./ vecnorm(sunRel,2,2));

x = (c.^2 + a.^2 - b.^2) ./ (2.*c);
y = sqrt(a.^2 - x.^2);
% y = real(sqrt(complex(a.^2 - x.^2))); % avoid being complex number

% occulted area
A = a.^2 .* acos(x./a) + b.^2 .* acos((c-x)./b) - c .* y; % avoid being complex number

nu = (a + b <= c) .* 1 ... % no occultation
    + (abs(a-b) < c) .* (c < a + b) .* (1 - A ./ pi ./ a.^2)... % partial
    + (c < b - a) .* (b > a) .* 0 ... % total occulation
    + (c < a - b) .* (a > b) .* (1 - b.^2 ./ a.^2); % total but maximum
end

%[appendix]{"version":"1.0"}
%---
