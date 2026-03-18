%[text] # convert day of year to Gregorian calendar date
%[text] ## inputs
%[text] `year`: year
%[text] `doy`: fractional day of year (1-based, 1.0 .. 366.x)
%[text] ## outputs
%[text] `month, day, hour, minute, second`
%[text] ## revisions
%[text] 20260317  y.yoshimura
%[text] See also gc2jd, jd2gc.
function [month, day, hour, minute, second] = doy2gc(year, doy)

jd = gc2jd(year, 1, 1, 0, 0, 0) + doy - 1;
[~, month, day, hour, minute, second] = jd2gc(jd);

end
