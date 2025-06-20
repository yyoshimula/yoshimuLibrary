%[text] # converting km to astronomical unit
%[text] `km`: distance to be converted, km
%[text] `const`: constant parameters for orbit propagation
%[text] `AU:` distance converted, AU
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20150101  y.yoshimura, y.yoshimula@gmail.com
%[text] See also au2km, orbitConst.
function AU = km2au(km, const)

AU = km ./ const.AU;

end

%[appendix]{"version":"1.0"}
%---
