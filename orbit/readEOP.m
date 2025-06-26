%[text] # read EOP data
%[text] ## note
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.213
%[text] ## revisions
%[text] 20230608  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function eopDataAll = readEOP(fName)

%[text] ## read EOP dat
eopDataAll = importdata(fName);


end



%[appendix]{"version":"1.0"}
%---
