%[text] # converting symbolic equations to latex script for Scrapbox
%[text] Scrapbox用にsymbolic equationを変換
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20220612  y.yoshimura, y.yoshimula@gmail.com
%[text] See also .
function out = sb(symEq)

out{:} = ['[$ ', latex(symEq), ']'];
fprintf('%s\n', out{:});

end

%[appendix]{"version":"1.0"}
%---
