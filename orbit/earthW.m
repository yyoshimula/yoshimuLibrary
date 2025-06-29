%[text] # Earth angular rate using LOD (length of day)
%[text] lod: length of day, nx1 vector, s
%[text] w: angular rate norm, nx1 vector, rad/s
%[text] ## note
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.222
%[text] ## revisions
%[text] 20230608  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function w = earthW(lod)
arguments
    lod (:,1) {mustBeNumeric}
end

w = 7.292115146706979 * 10^(-5) .* (1 - lod ./ 86400);

end



%[appendix]{"version":"1.0"}
%---
