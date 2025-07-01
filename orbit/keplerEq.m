%[text] # solving Kepler's equations
%[text] calculate eccentric anomaly from mean anomaly and eccentricity 
%[text] `M`: mean anomaly, rad, nx1
%[text] `e`: eccentricity nx1
%[text] TOL: tolerance
%[text] `E`: eccentric anomaly
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Curtis, Howard D 2013 Orbital mechanics for engineering students, p148 Algorithm 3.1 
%[text] ## revisions
%[text] 20211027  y.yoshimura, y.yoshimula@gmail.com
%[text] See also trueAnomaly.
function E = keplerEq(M, e, TOL)
% arguments
%     M (:,1)
%     e (:,1)
%     TOL = 1e-8
% end
if nargin < 3
    TOL = 1e-8;
end
residual = 1;
%[text] ## Kepler's equation
%[text] $M=E-e\\sin{E}$
%[text] Mによって初期推定値を変更
E = (M + 0.5 * e) .* (M < pi) + M - 0.5 * e .* (M > pi); % initial estimate

while(abs(residual) > TOL)
    fE = E - e * sin(E) - M;
    dfE = 1 - e * cos(E);

    residual = fE ./ dfE;
    E = E - residual;
end

%[text] `fsolve`を使った方法↓
% fun = @(x) x - e * sin(x) - M;
% E = fsolve(fun, E);

%[text] $E\\in \[0,2\\pi\]$
E = mod(E, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
