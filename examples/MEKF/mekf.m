function [qOut, xOut, pOut] = mekf(dt, xEst, P, qEst, qObs, wObs, Qnu, Qu, R)
% quaternionを観測量とする場合

dAlp = xEst(1:3);
dBet = xEst(4:6);

%% propagation step
Q = [(Qnu^2*dt+1/3*Qu^2.*dt^3).*eye(3), (1/2.*Qu^2.*dt^2).*eye(3)
    1/2.*Qu^2.*dt^2.*eye(3), Qu^2.*dt.*eye(3)];

% error state transition matrix, Phi
wEst = wObs - dBet;

wNorm = norm(wEst);
Phi11 = eye(3) - skew(wEst) .* (sin(wNorm .* dt))./wNorm ...
    + skew(wEst) * skew(wEst) * (1 - cos(wNorm * dt)) / wNorm^2;
Phi12 = skew(wEst) * (1 - cos(wNorm * dt)) / wNorm^2 - eye(3) .* dt ...
    - skew(wEst) * skew(wEst) * (wNorm * dt - sin(wNorm * dt)) / wNorm^3;
Phi21 = zeros(3,3);
Phi22 = eye(3);

Phi = [Phi11, Phi12
    Phi21, Phi22];

tmp = [-eye(3), zeros(3,3)
    zeros(3,3) eye(3)];

P = Phi * P * Phi' + tmp * Q * tmp';

% (global) quaternion propagation
co = cos(0.5 * norm(wEst) * dt);
si = sin(0.5 * norm(wEst) * dt);
nVec = wEst ./ norm(wEst);
qw = [nVec .* si,  co];
qEst = qMult(4,1, qw, qEst);

%% observation step
H = [eye(3), zeros(3,3)];

K = P * H' * (H * P * H' + R)^(-1);
P = (eye(6) - K * H) * P;

yTmp = qErr(4, qObs, qEst); % 4x1
yTmp = yTmp ./ norm(yTmp);
y = 2 .* yTmp(1:3) ./ yTmp(4);
y = y(:); % 3x1

tmp = K * (y -  dAlp');
xEst = tmp'; % 1x6

%% reset
qOut = qMult(4, 1, [xEst(1:3), 1], qEst);
qOut = qOut ./ norm(qOut);
xOut = [0, 0, 0, xEst(4:6)];
pOut = P;

end

