% solving Kepler's equations
% calculate eccentric anomaly from mean anomaly and eccentricity
% M: mean anomaly, rad, nx1
% e: eccentricity nx1
% E: eccentric anomaly
% note
% NA
% references
% Curtis, Howard D 2013 Orbital mechanics for engineering students, p148 Algorithm 3.1
% revisions
% 20211027  y.yoshimura, y.yoshimula@gmail.com
% See also trueAnomaly.

function E = keplerEq_(M, e)

TOL = 1e-8; % tolerance
residual = 1;

% Kepler's equation

% Mによって初期推定値を変更
E = (M + 0.5 * e) * (M < pi) + M - 0.5 * e * (M > pi); % initial estimate

while(abs(residual) > TOL)
    fE = E - e * sin(E) - M;
    dfE = 1 - e * cos(E);

    residual = fE / dfE;
    E = E - residual;
end

% fsolveを使った方法↓
% fun = @(x) x - e * sin(x) - M;
% E = fsolve(fun, E);

E = mod(E, 2*pi);

end
