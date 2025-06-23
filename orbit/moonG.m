%[text] # Moon gravitational force
%[text] ## inputs
%[text] `jd`: Julian day, day, 
%[text] `rVec`: satellite position at inertial frame, km, nx3 vector
%[text] `const`: orbital constants
%[text] `earthVSOP`, earth VSOP87 constants
%[text] `ELP`: ELP coefficients (not necessary if SPICE is used)
%[text] `options`
%[text]  `SPICE`: use SPICE, `on` or `off`
%[text] ## outputs
%[text] `aMoon`: Moon gravitational force, km/s^2, nx3 vector
%[text] `moonIJK`: moon position at inertial frame, km, nx3 vector
%[text] ## note
%[text] ELP is not necessary if SPICE is used
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, sunG.
function [aMoon, moonIJK] = moonG(jd, rVec, const, ELP, options)
% arguments
%     jd (:,1) {mustBeNumeric}
%     rVec (:,3) {mustBeNumeric}
%     const
%     ELP
%     options.SPICE char {mustBeMember(options.SPICE,{'on', 'off'})} = 'off'
% end
if nargin < 5
    options.SPICE = 'off';
end

%[text] ### default: SPICE is not used
if strcmp(options.SPICE, 'off')
    moonIJK = moon(jd, const, ELP);

else % when SPIC is used
    %     Julian day to ephemeris time (et)
    et = cspice_unitim(jd, 'JDTDT', 'ET');

    [moonIJK, ~] = cspice_spkpos('MOON', et, 'J2000', 'LT+S', 'EARTH');
    moonIJK = moonIJK'; % nx3
end

%[text] ### relative position vector from satellite to moon at inertial frame
rSC2m = moonIJK - rVec; % nx3

%[text] ### acceleration at inertial frame
aMoon = const.GM .* (rSC2m ./ vecnorm(rSC2m,2,2).^3 - moonIJK ./ vecnorm(moonIJK,2,2).^3);

end

%[appendix]{"version":"1.0"}
%---
