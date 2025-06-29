%[text] # DCM from PEF to ITRF
%[text] `xp, yp`: Earth orientation parameters, rad
%[text] R: direction cosine matrix, 3x3
%[text] ## note
%[text] this function is the rotation from PEF to ITRF (opposed to the Vallado's textbook)
%[text] ValladoはITRFからPEFへの変換にしているが，この関数はPEF to ITRF．
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p223. Eq. (3-78)
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also orbitConst, precession.
function R = pef2itrf(xp, yp)
% arguments
%     xp (1,1) {mustBeNumeric}
%     yp (1,1) {mustBeNumeric}    
% end


R = [1 0 xp
    0 1 -yp
    -xp yp 1];


end

%[appendix]{"version":"1.0"}
%---
