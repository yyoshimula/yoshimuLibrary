%[text] # calculating insertion date of leapseconds
%[text] `leapJD`:  （合計の）うるう秒が有効になるJulian dayと（合計の）うるう秒
%[text] ## note
%[text] `leapS.mを実行し，その結果をdAT.mの引数leapJDに使う．`
%[text] `うるう秒は挿入日の末尾（23:59:60 UTC）に挿入されるため，新しいTAI-UTCはその翌日0時（UTC）から有効になる．leapJD(:,1)は有効開始日のJDである．`
%[text] ## references
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition.
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] 20260706  fixed effective date (day after insertion), y.yoshimura
%[text] See also dAT
function leapJD = leapS()

% time table
t = leapseconds;

y = year(t.Date); % year
m = month(t.Date); % month
d = day(t.Date); % day

%[text] ## leapsecond effective date (JD) and leapseconds
% the leap second is inserted at the end of t.Date (23:59:60 UTC), so the
% new TAI-UTC takes effect from 00:00 UTC on the following day (+1 day)
% leapJD = [gc2jd(y, m, d, 0, 0, 0) + 1 time2num(t.CumulativeAdjustment)]; % if Predictive% Maintenance Toolbox exists
leapJD = [gc2jd(y, m, d, 0, 0, 0) + 1 seconds(t.CumulativeAdjustment)];

end

%[appendix]{"version":"1.0"}
%---
