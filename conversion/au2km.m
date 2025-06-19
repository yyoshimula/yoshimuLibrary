%[text] # converting astronomical unit (AU) to km
%[text] `AU:` distance to be converted, AU
%[text] `const`: constant parameters for orbit propagation
%[text] `km`: distance converted, km
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also km2AU, orbitConst.
function km = au2km(AU, const)

km = AU .* const.AU;

end

%[appendix]{"version":"1.0"}
%---
