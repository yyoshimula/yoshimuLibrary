%[text] # read ELP coefficients for Moon
%[text] ## note
%[text] ## references 
%[text] 
%[text] ## revisions
%[text] 20250625  y.yoshimura, y.yoshimula@gmail.com
%[text] See also moonG
function ELP = readELP(fNames)

%[text] ## read EOP data
if nargin < 1 % default files
    fNames(1) = {'ELPcoeffA.csv'};
    fNames(2) = {'ELPcoeffB.csv'};
end

ELP.a = readmatrix(fNames{1});
ELP.b = readmatrix(fNames{2});

end




%[appendix]{"version":"1.0"}
%---
