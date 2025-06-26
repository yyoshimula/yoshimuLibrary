%[text] # reading EGM2008 coefficients
%[text] EGM2008の係数読み込み
%[text] `fName`: coefficient file name, string
%[text] `deg`: degree to be read, scalar, int
%[text] `varagin`: if 'normalized' is set, outputs are normalized coefficients
%[text] `Clm`: coefficients, C, (deg+1) x (deg+1) matrix
%[text] `Slm`: coefficients, S, (deg+1) x (deg+1) matrix
%[text] ## note
%[text] 引数で'normalized'を指定したときだけnormalized coefficientsを出力
%[text] egm2008.mlxではunnormalized coefficientsを使用する
%[text] ## references 
%[text] David A. Vallado, "Fundamentals of Astrodynamics and Applications, 4th edition, pp.538–550.
%[text] Montenbruck  Oliver  & Eberhard Gill, Satellite Orbits. Springer Science & Business Media  2012., pp.61–67
%[text] ## revisions
%[text] 20210502  y.yoshimura
%[text] See also egm2008.
function EGM = readEGM2008(EGM, deg, normalized)
%[text] 引数で `normalized`を指定したときだけnormalized coefficientsを出力
if nargin < 3
    normalized = 0; % unnormalized coefficients
end

% file name (can be changed if another file is used) 
fName = 'EGM2008_to2190_TideFree.txt';

Cnm = zeros(deg+1, deg+1); % matlabのindexは1始まりなので+1
Snm = zeros(deg+1, deg+1);

fID = fopen(fName, 'r');
%[text] ## reading
while true
    tmp = sscanf(fgetl(fID), '%d%d%lfD%d%lfD%d%*lf%*lf'); %1行ずつ読み取り
    n = tmp(1);
    m = tmp(2);
    if n > deg
        break;
    end
    Cnm(n+1,m+1) = tmp(3)*10^(tmp(4));
    Snm(n+1,m+1) = tmp(5)*10^(tmp(6));

    if normalized == 0 % trnasform to unnormalized coefficients
        if m == 0
            k = 1;
        else
            k = 2;
        end
        Pi = sqrt(factorial(n+m)/(factorial(n-m) * k * (2*n+1)));
        Cnm(n+1,m+1) = Cnm(n+1,m+1) / Pi;
        Snm(n+1,m+1) = Snm(n+1,m+1) / Pi;
    end
end

fclose(fID);

EGM.Cnm = Cnm;
EGM.Snm = Snm;

end

%[appendix]{"version":"1.0"}
%---
