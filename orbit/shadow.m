%[text] # Earth shadow function
%[text] `satI`: satellite position@inertial frame, arbitrary unit, nx3
%[text] `sunI`: sun position@inertial frame, nx3
%[text] `rS`: Sun radius, scalar
%[text] `rE`: Earth radius, scalar
%[text] nu: 1: sunlit, 0: umbra (eclipse)
%[text] ## note
%[text] all variables (`satI, sunI, rS, rE`) must have the same unit
%[text] 変数の単位は統一すること．（統一していればmでもkmでもok）
%[text] ## references 
%[text] Satellite Orbits, Montenbruck, Gill, p.81
%[text] ## revisions
%[text] 20231128 y.yoshimura: 引数を変更
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also vsopConst.
function nu = shadow(satI, sunI, rS, rE)

sunRel = sunI - satI;

a = asin(rS ./ vecnorm(sunRel,2,2));  % Sun angular radius
b = asin(rE ./ vecnorm(satI,2,2));    % Earth angular radius
c = acos( dot(-satI, sunRel,2) ./ (vecnorm(satI,2,2).*vecnorm(sunRel,2,2)) ); % separation

n = size(satI,1);
nu = ones(n,1);

% total occultation (Sun fully covered): c <= b - a
idx_total = (c <= (b - a));
nu(idx_total) = 0;

% annular-like (Earth fully inside Sun disk): c <= a - b
idx_annular = (c <= (a - b)) & ~idx_total;
nu(idx_annular) = 1 - (b(idx_annular).^2 ./ a(idx_annular).^2);

% partial overlap: |a-b| < c < a+b
idx_part = (abs(a-b) < c) & (c < (a+b));
if any(idx_part)
    ai = a(idx_part); 
    bi = b(idx_part); 
    ci = c(idx_part);

    u = (ci.^2 + ai.^2 - bi.^2) ./ (2*ci);
    v = (ci.^2 + bi.^2 - ai.^2) ./ (2*ci);

    % clip for numerical safety
    t1 = max(-1, min(1, u./ai));
    t2 = max(-1, min(1, v./bi));

    w  = max(0, ai.^2 - u.^2);

    A = ai.^2 .* acos(t1) + bi.^2 .* acos(t2) - ci .* sqrt(w);
    nu(idx_part) = 1 - A ./ (pi*ai.^2);
end

% no overlap already nu=1: c >= a+b

end


%[appendix]{"version":"1.0"}
%---
