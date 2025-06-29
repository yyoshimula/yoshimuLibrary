%[text] # quaternion from ITRF to GCRF (CIO approcach of　IAU-2006/2000 reduction)
%[text] ## inputs
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$\]
%[text] R, rotation matrix, 3x3 matrix 
%[text] ## note
%[text] itrf2gcrfのwrapper的な関数．使いやすいように引数を変更している．
%[text] 位置座標の変換のみ．速度を変換する際は地球の自転速度を考慮した計算が必要
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.220
%[text] ## revisions
%[text] 20230612  y.yoshimura, y.yoshimula@gmail.com
%[text] See also itrf2gcrf.
function q = qITRF2gcrf(scalar, jd)
% arguments
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     jd
% end

% for time calculation
iau06 = readIAU06;
leapJD = leapS; % load database
fname = "EOP_20_C04_one_file_1962-now.txt";
eopDataAll = readEOP(fname); % earth oriented parameters (EOP)

[yy, mm, dd, hh, m, s] = jd2gc(jd);
utc = datetime([yy, mm, dd, hh, m, s]); % UTC (datetime format)

for i = 1:length(jd)
    dcmITRF2gcrf = itrf2gcrf(utc(i), leapJD, iau06, eopDataAll);
    q(i,:) = dcm2q(scalar, dcmITRF2gcrf);
end


%[appendix]{"version":"1.0"}
%---
