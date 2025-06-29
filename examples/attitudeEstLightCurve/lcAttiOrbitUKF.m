%[text] # photometric attittude estimtion using UKF
%[text] state vector
%[text] x:= \[positiona, error GRP\]
function [xEstOut, qGlobalOut, PcovOut] = lcAttiOrbitUKF(jd, dt, xEst, Pcov, y, aziObs, eleObs, ...
    ukfPara, Qest, Rest, fGRP, aGRP, qGlobal, w, v, obsECI, sat, const, earthVSOP)

n_ = length(xEst);

[x0, xx] = ukfSigma(ukfPara.lam, Pcov, xEst);

% error quaternion sigma points
qe0Sigma = grp2q(4, fGRP, aGRP, x0(4:6));
qeSigma = grp2q(4, fGRP, aGRP, xx(:,4:6)); % 2n x 4

% global quaternion sigma points
q0Sigma = qGlobal;
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % 2nx4
%[text] ## propagation step
%[text] ### attitude
yTmp = qPropMat(4, dt, w) * q0Sigma';  % closed form

q0Sigma = yTmp';
qe0Sigma = qErr(4, q0Sigma, q0Sigma);
x0(4:6) = q2GRP(4, fGRP, aGRP, qe0Sigma);

for j = 1:2*n_
    yTmp = qPropMat(4,dt, w) * qSigma(j,:)';
    qSigma(j,:) = yTmp';

    % add process noise
    Qq = grp2q(4, fGRP, aGRP, diag(Qest)');
    qSigma(j,:) = qMult(4, 1, Qq, qSigma(j,:));

    % quaternion to error GRP
    qeSigma(j,:) = qErr(4, qSigma(j,:), q0Sigma);
    xx(j,4:6) = q2GRP(4, fGRP, aGRP, qeSigma(j,:)); % error GRP sigma point
end
%[text] ### orbit
[~, tmp] = ode113(@(t,x)eomOrbit(t,x,const), 0:0.1:5, [x0(1:3), v]);
x0(1:3) = tmp(end,1:3);
for j = 1:2*n_
    [~, tmp] = ode113(@(t,x) eomOrbit(t,x,const),  0:0.1:5, [xx(j,1:3), v]);
    xx(j,1:3) = tmp(end,1:3);
end
%[text] ### mean and covariance
xEst = ukfPara.w0m .* x0 + ukfPara.wim .* sum(xx,1); % predicted state
Pcov = ukfCov(xEst, x0, xx, ukfPara.w0c, ukfPara.wic, 0);
%[text] ### recalculate sigma points
[x0, xx] = ukfSigma(ukfPara.lam, Pcov, xEst);

% error quaternion sigma points
qe0Sigma = grp2q(4, fGRP, aGRP, x0(4:6));
qeSigma = grp2q(4, fGRP, aGRP, xx(:,4:6)); % 2n x 4

% global quaternion sigma points
q0Sigma = qGlobal;
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % 2nx4
%[text] ## predicted output sigma points
[lon, lat, dist] = sun(jd, const, earthVSOP);
[sunPos(1), sunPos(2), sunPos(3)] = sph2cart(lon, lat, au2km(dist, const)); % sun position @inertial frame
sunPos = sunPos .* 10^3; % m
obsPos = obsECI; % observer position @inertial frame, m

if (y ~= Inf) % observation updateする
    nu = 1; % earth shadowing, とりあえず
    satPos = x0(1:3); % sat position @inertial frame, m
    dRel = norm(obsPos - satPos); % m
    [~, fObs] = lc(sat, 4, q0Sigma, satPos, obsPos, sunPos, 1, 'AS');
    ye0 = fObs;
    azi0 = atan2(satPos(2), satPos(1));
    ele0 = asin(satPos(3) / norm(satPos));

    for j = 1:2*n_
        satPos = xx(j,1:3); % sat position @inertial frame, m
        dRel = norm(obsPos - satPos); % m
        [~, fObs] = lc(sat, 4, qSigma(j,:), satPos, obsPos, sunPos, 1, 'AS');
        yez(j,1) = fObs;
        azi(j,1) = atan2(satPos(2), satPos(1));
        ele(j,1) = asin(satPos(3) / norm(satPos));
    end

    ye0 = mag(ye0, dRel);
    yez = mag(yez, dRel);

    ye0 = [ye0, azi0, ele0];
    yez = [yez, azi, ele];

    % predicted mean output
    ye = ukfPara.w0m .* ye0 + ukfPara.wim .* sum(yez, 1);
    if(isnan(ye))
        disp('no observation update');
        xEst = xEst;
        p = diag(Pcov)';
    else %observation update
        % Calculate correlation
        [Pyy, ~, K] = ukfCorrGain(xEst, x0, xx, ye, ye0, yez, ukfPara.w0c, ukfPara.wic, Rest);

        % Update
        Pcov = Pcov - K * Pyy * K';
        % m_appを観測値とする場合
        xEst = xEst + (K * ([y,aziObs,eleObs] - ye)')';
    end

else
    disp('no observation update')
end
%[text] ## Global quaternion and reset error GRP
% global quaternion
qeTmp = GRP2q(4, fGRP, aGRP, xEst(1,4:6));
qGlobalOut = qMult(4, 1, qeTmp, q0Sigma);

PcovOut = Pcov;

xEstOut = xEst;
% error GRP reset
xEstOut(4:6) = [0 0 0];
end

%[appendix]{"version":"1.0"}
%---
