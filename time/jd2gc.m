%[text] # calculate the Gregorian calendar date from the Julian date
%[text] `jd`: julian day, day
%[text] `year`
%[text] `month`
%[text] `day`
%[text] `hour`
%[text] `minute`
%[text] `second`
%[text] ## note
%[text] `checked with:` [`http://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/jd2date.cgi`](http://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/jd2date.cgi)
%[text] ## references
%[text] `NA`
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also jd2fyear.
function [year, month, day, hour, minute, second] = jd2gc(jd)

fday   = jd + 0.5 - floor(jd + 0.5);
jd     = floor(jd + 0.5);
in2    = floor(floor((jd - 4479.5) / 36524.25) * 0.75 + 0.5) - 37;
in2    = jd + in2;

year   = floor(in2 / 365.25) - 4712;
an2    = in2 - 59.25;
an2    = floor(an2 - floor(an2 / 365.25) * 365.25) + 0.5;
mon2   = floor(an2 / 30.6) + 2;

month  = mon2 - floor(mon2 / 12) * 12 + 1;
day    = floor(an2 - floor(an2 / 30.6) * 30.6 + 1);
hour   = floor(fday * 24);
minute = floor((fday * 24 - floor(fday * 24)) * 60);
second = fday * 24.0 * 60.0 * 60.0 - floor(fday * 24) * 60.0 * 60.0...
    - floor((fday * 24 - floor(fday * 24)) * 60) * 60.0;


end

%[appendix]{"version":"1.0"}
%---
