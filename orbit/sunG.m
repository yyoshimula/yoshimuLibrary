%[text] # Sun's gravitational force
%[text] ## input 
%[text] `jd`: Julian day, day, 
%[text] `rVec`: satellite position at inertial frame, km, nx3 vector
%[text] `const`: orbital constants
%[text] `earthVSOP`, earth VSOP87 constants
%[text] `options`
%[text]  `SPICE`: use SPICE, `on` or `off`
%[text] ## `output`
%[text] `aSun`: Sun's gravitational force, km/s^2, nx3 vector
%[text] `sunIJK`: Sun's position at inertial frame, km, nx3 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, vsopConst, precession.
function [aSun, sunI] = sunG(jd, rVec, const, earthVSOP, options)
% arguments
%     jd (:,1) {mustBeNumeric}
%     rVec (:,3) {mustBeNumeric}
%     const
%     earthVSOP
%     options.SPICE char {mustBeMember(options.SPICE,{'on', 'off'})} = 'off'
% end
if nargin < 5
    options.SPICE = 'off';
end

%[text] ### default: SPICE is not used
if strcmp(options.SPICE, 'off')
    sunI = sun(jd, const, earthVSOP);

else % when SPIC is used
    %     Julian day to ephemeris time (et)
    et = cspice_unitim(jd, 'JDTDT', 'ET');

    [sunI, ~] = cspice_spkpos('SUN', et, 'J2000', 'LT+S', 'EARTH');
    sunI = sunI';
end

% relative position vector from satellite to sun at inertial frame
rSC2s = sunI - rVec; % nx3

% at inertial frame
aSun = const.GS .* (rSC2s ./ vecnorm(rSC2s,2,2).^3 - sunI ./ vecnorm(sunI,2,2).^3);

end

%[appendix]{"version":"1.0"}
%---
