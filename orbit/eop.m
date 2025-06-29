%[text] # get earth orientation parameters (EOP) for IAU-2000 reduction
%[text] EOPを抽出
%[text] `mjd:` modified julian day
%[text] `xp, yp:` $x\_p, y\_p\n\n$ for DCM of polar motion 
%[text] `dUT1:` $\\rm \\Delta UT1$
%[text] `dX, dY: for DCM of polar motion`
%[text] `lod:` length of day
%[text] ## note
%[text] IERSから抽出してコメント部分を消したデータが必要
%[text] ここではサンプルとしてEOP\_20\_C04\_one\_file\_1962-now.txtを使用．
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition.
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function output = eop(year, month, day, eopDataAll)
% arguments
%     year (:,1) {mustBeNumeric}
%     month (:,1) {mustBeNumeric}
%     day (:,1) {mustBeNumeric}
%     eopDataAll
% end

% julian day to modified julian day
jd = gc2jd(year, month, day, 0, 0, 0);
mjd = jd2mjd(jd);

% find index
ind = find(eopDataAll.data(:,5) == mjd);

output.mjd = eopDataAll.data(ind,5);
output.xp = arcs2rad(eopDataAll.data(ind,6));
output.yp = arcs2rad(eopDataAll.data(ind,7));
output.dUT1 = eopDataAll.data(ind,8);
output.dX = arcs2rad(eopDataAll.data(ind,9));
output.dY = arcs2rad(eopDataAll.data(ind,10));
output.lod = eopDataAll.data(ind,13);

end


%[appendix]{"version":"1.0"}
%---
