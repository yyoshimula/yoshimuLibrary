%[text] # Gauss method of orbit determination
%[text] `t:` time at observations, 3x1 vector
%[text] `azi, ele:` azimuth and elevation angles of the object at topocentric equatorial frame, 3x2 matrix
%[text] `obsECI:` observer position vector at inertial frame, 3x3 matrix
%[text] `rRange`: initial guess range for the 8th order equations
%[text] `mu`: Earth's gravitational constant, scalar
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Curtis, Howard D. Orbital Mechanics for Engineering Students. Butterworth-Heinemann, 2013, p268.
%[text] ## revisions
%[text] 20221110  y.yoshimura, y.yoshimula@gmail.com, y.yoshimura.a64@m.kyushu-u.ac.jp
%[text] See also gibbs.
function [r2,v2] = gauss(t, aziele, obsECI, rRange, mu)

TOL = 1e-8; % tolerance for 8th order equation
residual = 1;
%[text] ## time intervals
t1 = t(1);
t2 = t(2);
t3 = t(3);

tau1 = t1 - t2;
tau3 = t3 - t2;
tau = tau3 - tau1;

%[text] ## unit vectors
[tmp1, tmp2, tmp3] = sph2cart(aziele(:,1), aziele(:,2), 1);
tmp = [tmp1, tmp2, tmp3];

rho1 = tmp(1,:);
rho2 = tmp(2,:);
rho3 = tmp(3,:);

p1 = cross(rho2, rho3);
p2 = cross(rho1, rho3);
p3 = cross(rho1, rho2);

%[text] ## calculate D
R1 = obsECI(1,:);
R2 = obsECI(2,:);
R3 = obsECI(3,:);

D0 = rho1 * p1';

D11 = R1 * p1';
D12 = R1 * p2';
D13 = R1 * p3';

D21 = R2 * p1';
D22 = R2 * p2';
D23 = R2 * p3';

D31 = R3 * p1';
D32 = R3 * p2';
D33 = R3 * p3';

%[text] ## calculate A and B (in Eq.(5.112a))
A = (-D12 * tau3 / tau + D22 + D32 * tau1 / tau) / D0;
B = (D12 * (tau3^2 - tau^2) * tau3 / tau + D32 * (tau^2 - tau1^2) * tau1 / tau) / 6 / D0;

%[text] ## calculate E, a, b, and c
E = R2 * rho2';

a = -(A^2 + 2 * A * E + norm(R2)^2);
b = -2 * mu * B * (A + E);
c = -mu^2 * B^2;

figure
fx = rRange.^8 + a .* rRange.^6 + b .* rRange.^3 + c;
plot(rRange, fx)

%[text] ### figureを表示して，そこからinitial guessを設定
% ライン上をマウスで選択
waitforbuttonpress; % 入力待ち状態

% マウスの現在値の取得
Cp = get(gca,'CurrentPoint');
Xp = Cp(2,1);  % X座標
Yp = Cp(2,2);  % Y座標

% マウスで選択した位置と最も近いデータのインデックスを計算
[~, Ip] = min((rRange-Xp).^2+(fx-Yp).^2);
% マウスで選択した位置と最も近いデータ値を取得．それをinitial guessにする．
Xc = rRange(Ip);
% Yc = fx(Ip);
xIni = Xc;
%[text] ## solve 8th order equation
x = xIni;
while(abs(residual) > TOL)
    fx = x^8 + a * x^6 + b * x^3 + c;
    dfx = 8 * x^7 + 6 * a * x^5 + 3 * b * x^2;

    residual = fx / dfx;
    x = x - residual;
end
r2Norm = x;
%[text] ## slant ranges
rho1Norm = (6 * (D31 * tau1 / tau3 + D21 * tau / tau3) * r2Norm^3 + mu * D31 * (tau^2 - tau1^2) * tau1 / tau3) / (6 * r2Norm^3 + mu * (tau^2 - tau3^2)) - D11;
rho1Norm = rho1Norm / D0;

rho2Norm = A + mu * B / r2Norm^3;

rho3Norm = (6 * (D13 * tau3 / tau1 - D23 * tau / tau1) * r2Norm^3 + mu * D13 * (tau^2 - tau3^2) * tau3 / tau1) / (6 * r2Norm^3 + mu * (tau^2 - tau1^2)) - D33;
rho3Norm = rho3Norm / D0;

%[text] ## position vectors at inertial frame
r1 = R1 + rho1Norm .* rho1;
r2 = R2 + rho2Norm .* rho2;
r3 = R3 + rho3Norm .* rho3;

%[text] ## Lagrange coefficients
f1 = 1 - 0.5 * mu / r2Norm^3 * tau1^2;
f3 = 1 - 0.5 * mu / r2Norm^3 * tau3^2;

g1 = tau1 - mu / 6 / r2Norm^3 * tau1^3;
g3 = tau3 - mu / 6 / r2Norm^3 * tau3^3; 

%[text] ## velocity vectors
v2 = 1 / (f1 * g3 - f3 * g1) * (-f3 .* r1 + f1 .* r3);

end

%[appendix]{"version":"1.0"}
%---
