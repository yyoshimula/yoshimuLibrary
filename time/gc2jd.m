%[text] # converting Gregorian calendar date to Julian days
%[text] `year`: year, nx1 vector
%[text] `month: month,` nx1 vector
%[text] `hour`: hour`,` nx1 vector
%[text] `minute`: minute`,` nx1 vector
%[text] `second`: second`,` nx1 vector
%[text] `jd`: Julian day`,` nx1 vector
%[text] ## note
%[text] NA
%[text] ## references 
%[text] confirmed with [http://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/date2jd.cgi](http://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/date2jd.cgi)
%[text] Meeus, J., Astronomical Algorithms, 1998., p61, Eq.(7.1)
%[text] ## revisions
%[text] 20160630  y.yoshimura, y.yoshimula@gmail.com
%[text] See also jd2gc.
function jd = gc2jd(year, month, day, hour, minute, second)
% arguments
%     year (:,1) {mustBeNumeric}
%     month (:,1) {mustBeNumeric}
%     day (:,1) {mustBeNumeric}  
%     hour (:,1) {mustBeNumeric}
%     minute (:,1) {mustBeNumeric}
%     second (:,1) {mustBeNumeric}  
% end

M_ = (month > 2) .* month + (month <= 2) .* (month + 12);
Y_ = (month > 2) .* year  + (month <= 2) .* (year  -  1 );
 
A_ = floor(Y_ ./ 100);
B_ = 2 - A_ + floor(A_ ./ 4);

jd = floor(365.25 .* (Y_ + 4716)) + floor(30.6001 .* (M_ + 1))...
    + day + B_ - 1524.5 + hour ./ 24 + minute ./ 1440 + second ./ 86400;

end

%[appendix]{"version":"1.0"}
%---
