%[text] # read EOP data
%[text] ## note
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.213
%[text] ## revisions
%[text] 20230608  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function EOP = readEOP(fName)

%[text] ## read EOP dat
if nargin < 1
    fName = "EOP_20_C04_one_file_1962-now.txt";
end

EOP.dataAll = importdata(fName);

EOP.iau06 = readIAU06;
EOP.leapJD = leapS; % load database

end



%[appendix]{"version":"1.0"}
%---
