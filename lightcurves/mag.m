%[text] # calculating apparent magnitude of object with distance d \[m\]
%[text] `fobs`: total observed reflectance flux
%[text] `d`: distance between object and observer, m
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20251208  y.yoshimura, y.yoshimula@gmail.com
function appMag = mag(fobs, d)

msun = -26.7; % apparent magnitude of Sun
appMag = msun - 2.5 .* log10(fobs ./ d.^2);


end

%[appendix]{"version":"1.0"}
%---
