%[text] # rotation matrix from ITRF to GCRF (CIO approcach of　IAU-2006/2000 reduction)
%[text] `jd`:  julian day
%[text] `leapJD`: obtained from `leapS.mlx`
%[text] ## note
%[text] 位置座標の変換のみ．速度を変換する際は地球の自転速度を考慮した計算が必要
%[text] ## requirement
%[text] `leapS.mlx`を実行し，その出力`leapJD`を用いる 
%[text] `readIAU06.mlx`を実行し，その出力`iau06`を用いる （`iau06x.dat, iau06y.dat, iau06z.datを使用）`
%[text] `readEOP.mlx`を実行し，その出力`eopDataAll`を用いる 
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.220
%[text] ## revisions
%[text] 20230612  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gcrf2itrf, utc2tt.
function dcm = itrf2gcrf(jd, leapJD, iau06, eopDataAll)
% arguments
%     jd 
%     leapJD
%     iau06
%     eopDataAll
% end
 
%[text] ## time and EOP
[yy, mm, dd, hh, m, s] = jd2gc(jd);
utc = datetime([yy mm dd hh m s]);

deltaAT = dAT(jd, leapJD);
eopData = eop(yy, mm, dd, eopDataAll);

UT1 = utc + seconds(eopData.dUT1);

jdUT1 = gc2jd(year(UT1), month(UT1), day(UT1), hour(UT1), minute(UT1), second(UT1));

jdTT = utc2tt(jd, deltaAT);
tTT = jd2jdT(jdTT); % Julian century of TT

%[text] ## DCM (from ITRF to TIRS)
sPrime = wobble(tTT);
itrf2TIRS = dcm1axis(3, -sPrime) * dcm1axis(2, eopData.xp) * dcm1axis(1, eopData.yp);

%[text] #### TIRS
% rTIRS = itrf2TIRS * rITRF;
% vTIRS = itrf2TIRS * vITRF;

%[text] ## earth rotation angle
thetaERA = era(jdUT1);
%[text] ## DCM (from TIRS to CIRS)
tirs2CIRS = dcm1axis(3, -thetaERA);

%[text] memo
%[text] 速度ベクトルを非回転座標系へ変換する際は注意
% we = [0; 0; earthW(eopData.lod)];
% rCIRS = tirs2CIRS * rTIRS;
% vCIRS = tirs2CIRS * (vTIRS + cross(we,rTIRS))

%[text] ## precession and nutation (DCM from CIRS to GCRF)
cirs2GCRF = precessionNutation(jdTT, iau06, eopData.dX, eopData.dY);

dcm = cirs2GCRF * tirs2CIRS * itrf2TIRS;

end



%[appendix]{"version":"1.0"}
%---
