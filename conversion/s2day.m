%[text] # transform second to days 
%[text] ## input
%[text] s, second
%[text] ## output
%[text] day, day
%[text] 20210215 y.yoshimura 
function day = s2day(s)

day = s ./ 24 ./ 60 ./ 60;

end

%[appendix]{"version":"1.0"}
%---
