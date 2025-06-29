%[text] # annual wobble for polar motion of the Earth
%[text] # IAU-2006/2000, CIO-based
%[text] tTT: julian century
%[text] sPrime: rad
%[text] ## note
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.212
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function sPrime = wobble(tTT)
% arguments
%     tTT (:,1) {mustBeNumeric}
% end

sPrime = arcs2rad(-0.000047) * tTT;

end



%[appendix]{"version":"1.0"}
%---
