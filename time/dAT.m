%[text] # calculating ΔAT (= TAI - UTC)
%[text] dATを計算
%[text] ## note
%[text] `dAT.mlxを使う前にleapS.mlxを実行し，その結果を引数leapJDに使う．`
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition.
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function deltaAT = dAT(jd, leapJD)
% arguments
%     jd (:,1) {mustBeNumeric}
%     leapJD (:,2) {mustBeNumeric}
% end

%[text] 年代の新しい方から比較
i = length(leapJD);
while (i >= 1 && jd < leapJD(i,1))
    i = i - 1;
end

if i == 0 % before 1972-1-1
    deltaAT = 10;
else
    deltaAT = leapJD(i,2) + 10;
end

end

%[appendix]{"version":"1.0"}
%---
