%[text] # calculating ΔAT (= TAI - UTC)
%[text] dATを計算
%[text] ## note
%[text] `dAT.mを使う前にleapS.mを実行し，その結果を引数leapJDに使う．`
%[text] `jdはスカラでもベクトルでもよい（deltaATはjdと同じサイズで返る）．`
%[text] ## references
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition.
%[text] ## revisions
%[text] 20230605  y.yoshimura, y.yoshimula@gmail.com
%[text] 20260706  vectorized for jd input, y.yoshimura
%[text] See also leapS
function deltaAT = dAT(jd, leapJD)
% arguments
%     jd (:,1) {mustBeNumeric}
%     leapJD (:,2) {mustBeNumeric}
% end

%[text] 各jdに対して，有効開始JD（leapJD(:,1)）がjd以前のエントリ数を数える
n = sum(jd(:) >= leapJD(:,1).', 2);

% n = 0 (before the first leap second): TAI - UTC = 10 s (since 1972-1-1)
cumLeap = [0; leapJD(:,2)];
deltaAT = reshape(10 + cumLeap(n + 1), size(jd));

end

%[appendix]{"version":"1.0"}
%---
