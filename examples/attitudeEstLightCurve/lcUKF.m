function    [xEstOut, qOut, Pout] = lcUKF(jd, dt_, xEst, Pcov, y, ...
    ukfPara, Qest, Rest, fGRP, aGRP, qGlobal, w, r, obsECI, sat, const, earthVSOP)

n_ = length(xEst); % # of state variables

%% initial sigma points
[x0, xx] = ukfSigma(ukfPara.lambda, Pcov, xEst);

% error quaternion sigma points
qe0Sigma = grp2q(4, fGRP, aGRP, x0(1:3));
qeSigma = grp2q(4, fGRP, aGRP, xx(:,1:3)); % 2n x 4

% global quaternion sigma points
q0Sigma = qGlobal;
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % 2nx4

%% propagation step
% X0
tmp = qPropMat(4, dt_, w) * q0Sigma';  % discrete closed form
q0Sigma = tmp';
qe0Sigma = qErr(4, q0Sigma, q0Sigma);
x0(1:3) = q2GRP(4, fGRP, aGRP, qe0Sigma);

% X
for j = 1:2*n_
    tmp = qPropMat(4,dt_, w) * qSigma(j,:)';
    qSigma(j,:) = tmp';

    % add process noise
    % Qq = grp2q(4, fGRP, aGRP, diag(Qest)');
    % qSigma(j,:) = qMult(4, 1, Qq, qSigma(j,:));

    % quaternion to error GRP
    qeSigma(j,:) = qErr(4, qSigma(j,:), q0Sigma);
    xx(j,1:3) = q2GRP(4, fGRP, aGRP, qeSigma(j,:)); % error GRP sigma point
end
% mean and covariance
xEst = ukfPara.w0m .* x0 + ukfPara.wim .* sum(xx,1); % predicted state
Pcov = ukfCov(xEst, x0, xx, ukfPara.w0c, ukfPara.wic, 0 .* Qest);

% sigma point recalculation
[x0, xx] = ukfSigma(ukfPara.lambda, Pcov, xEst);

%% observation step
sunPos = sun(jd, const, earthVSOP); % km, sun position @inertial frame (J2000.0)
sunPos = sunPos .* 10^3; % m
satPos = r; % sat position @inertial frame, m
obsPos = obsECI; % observer position @inertial frame, m

if (y ~= Inf) % observation updateする
    nu = 1; % earth shadowing, とりあえず

    [mApp, ~] = lc(sat, 4, q0Sigma, satPos, obsPos, sunPos, nu, 'AS');
    ye0 = mApp;

    for j = 1:2*n_
        [mApp, ~] = lc(sat, 4, qSigma(j,:), satPos, obsPos, sunPos, nu, 'AS');
        yez(j,1) = mApp;
    end

    % predicted mean output
    ye = ukfPara.w0m .* ye0 + ukfPara.wim .* sum(yez, 1);

    if(isnan(ye) || isinf(ye))
        disp('no observation update');
        xEst = xEst;
        p = diag(Pcov)';

    else %observation update
        % Calculate correlation
        [Pyy, ~, K] = ukfCorrGain(xEst, x0, xx, ye, ye0, yez, ukfPara.w0c, ukfPara.wic, Rest);

        % Update
        Pcov = Pcov - K * Pyy * K';
        p = diag(Pcov)';
        
        xEst = xEst + (K * (y - ye))';
    end
else
    disp('no observation update')

end

%% Global quaternion and reset error GRP
% global quaternion
qeTmp = grp2q(4, fGRP, aGRP, xEst(1,1:3));
qOut = qMult(4, 1, qeTmp, q0Sigma);

Pout = Pcov;

% error GRP reset
xEstOut = zeros(1,3);
end
