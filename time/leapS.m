%[text] # calculating insertion date of leapseconds
%[text] `leapJD`:  うるう秒が挿入されたJulian dayと（合計の）うるう秒
%[text] ## note
%[text] `leapS.mを実行し，その結果をdAT.mの引数leapJDに使う．`
%[text] ## references
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition.
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also dAT
function leapJD = leapS()

% time table
t = leapseconds;

y = year(t.Date); % year
m = month(t.Date); % month
d = day(t.Date); % day

%[text] ## leapsecond date (JD) and leapseconds
% leapJD = [gc2jd(y, m, d, 0, 0, 0) time2num(t.CumulativeAdjustment)]; % if Predictive% Maintenance Toolbox exists
leapJD = [gc2jd(y, m, d, 0, 0, 0) seconds(t.CumulativeAdjustment)];

end

%[appendix]{"version":"1.0"}
%---
