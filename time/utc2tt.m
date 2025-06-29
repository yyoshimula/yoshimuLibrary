%[text] # UTC to TT (Terrestrial time)
%[text] ## inputs
%[text] `jdUTC: Julian day, day`
%[text] `deltaAT: dAT, s`
%[text] ## outputs
%[text] `jdTT`: julian date of TT, day
%[text] ## note
%[text] $TT = TAI +32.184\\\\\n=UTC+\\Delta AT +32.184 $
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. p.220
%[text] ## revisions
%[text] 20230612  y.yoshimura
%[text] See also orbitConst
function jdTT = utc2tt(jdUTC, deltaAT)
% arguments
%     jdUTC (:,1) {mustBeNumeric}
%     deltaAT (:,1) {mustBeNumeric}
% end

jdTAI = jdUTC + s2day(deltaAT);
jdTT = jdTAI + s2day(32.184);

end

%[appendix]{"version":"1.0"}
%---
