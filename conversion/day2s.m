%[text] # transform days to seconds
%[text] ## input
%[text] day, day
%[text] ## output
%[text] s, second
%[text] ## note
%[text] `s2day` の逆変換。日単位の時間差を秒に変換する。
%[text] See also s2day.
%[text] ## revisions
%[text] 20260707 y.yoshimura
function s = day2s(day)

s = day .* 24 .* 60 .* 60;

end

%[appendix]{"version":"1.0"}
%---
