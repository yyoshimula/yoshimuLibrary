%[text] # weighted circular mean of angles
%[text] 角度の重み付き円周平均を計算 (単純な算術平均ではなく、cos/sin の重み付き和を atan2 で合成)
%[text] ## inputs
%[text] `thetaArray`: angles, rad, vector
%[text] `w`: weights, vector (same length as thetaArray, normalized internally)
%[text] ## output
%[text] `thetaAve`: weighted mean angle, $\in (-\pi, \pi]$, rad, scalar
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
function thetaAve = meanAngle(thetaArray, w)
    thetaArray = thetaArray(:);
    w = w(:);
    w = w ./ sum(w);

    x = sum(w .* cos(thetaArray));
    y = sum(w .* sin(thetaArray));
    
    thetaAve = atan2(y, x);
end