%[text] # Earth rotation angle (IAU-2006/2000, CIO-based)
%[text] jdUT1: julian day of UT1, day
%[text] theta: earth rotation angle, rad
%[text] ## note
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.212
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function theta = era(jdUT1)
% arguments
%     jdUT1 (:,1) {mustBeNumeric}
% end

theta = 0.779057273264 + 1.00273781191135448 * (jdUT1 - 2451545.0);
theta = mod(theta .* 2.0 .* pi, 2*pi);

end



%[appendix]{"version":"1.0"}
%---
