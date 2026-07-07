%[text] # degree, arcminute, arcsecond to degree (DMS → deg)
%[text] 度・分・秒 (DMS) 表記の角度を度に変換する。
%[text] ## input
%[text] hour: degree part, deg
%[text] min: arcminute, arcmin
%[text] sec: arcsecond, arcsec
%[text] ## output
%[text] out: angle, deg
%[text] ## note
%[text] deg = hour + min/60 + sec/3600（RA 時系の ×15 換算ではない）。唯一の呼び出し元 obliquity.m の黄道傾斜角 23°26'21.448" と整合する。
%[text] See also obliquity.
%[text] ## revisions
%[text] 20260707 y.yoshimura, hour 項の誤った ×15 を除去し DMS→deg として自己整合化
%[text] 20210419 y.yoshimura
function out = hms2deg(hour, min, sec)

out = hour + (min .* 60 + sec) ./ 3600;


end

%[appendix]{"version":"1.0"}
%---
