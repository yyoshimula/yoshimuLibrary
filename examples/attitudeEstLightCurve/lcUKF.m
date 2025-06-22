function    [xEstOut, qOut, Pout] = lcUKF(jd, dt_, xEst, Pcov, y, ...
    ukfPara, Qest, Rest, fGRP, aGRP, qGlobal, w, r, obsECI, sat, const, earthVSOP)

n_ = length(xEst); % # of state variables

%% initial sigma points
xx = ukfSigma(ukfPara.lambda, Pcov, xEst);

% error quaternion sigma points
qeSigma = grp2q(4, fGRP, aGRP, xx(:,1:3)); % 2n+1 x 4

% global quaternion sigma points
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % 2n+1 x 4

%% propagation step
% X
for j = 1:(2*n_+1)
    tmp = qPropMat(4,dt_, w) * qSigma(j,:)';
    qSigma(j,:) = tmp';

    % add process noise
    % Qq = grp2q(4, fGRP, aGRP, diag(Qest)');
    % qSigma(j,:) = qMult(4, 1, Qq, qSigma(j,:));

    % quaternion to error GRP
    qeSigma(j,:) = qErr(4, qSigma(j,:), qSigma(1,:));
    xx(j,1:3) = q2grp(4, fGRP, aGRP, qeSigma(j,:)); % error GRP sigma point
end
% mean and covariance
xEst = sum(ukfPara.wm .* xx, 1); % predicted state
Pcov = ukfCov(xEst, xx, ukfPara.wc, 0 .* Qest);

% sigma point recalculation
xx = ukfSigma(ukfPara.lambda, Pcov, xEst);

%% observation step
sunPos = sun(jd, const, earthVSOP); % km, sun position @inertial frame (J2000.0)
sunPos = sunPos .* 10^3; % m
satPos = r; % sat position @inertial frame, m
obsPos = obsECI; % observer position @inertial frame, m

options.BRDF = "AS";
options.mex = "off";

if (y ~= Inf) % observation updateする
    nu = 1; % earth shadowing, とりあえず

    for j = 1:(2*n_+1)
        [mApp, ~] = lc(sat, 4, qSigma(j,:), satPos, obsPos, sunPos, nu, options);
        yez(j,1) = mApp;
    end

    % predicted mean output
    yEst = sum(ukfPara.wm .* yez, 1);

    if(isnan(yEst) || isinf(yEst))
        disp('no observation update');

    else %observation update
        % Calculate correlation
        [Pyy, ~, K] = ukfCorrGain(xEst, xx, yEst, yez, ukfPara.wc, Rest);

        % Update
        Pcov = Pcov - K * Pyy * K';
        p = diag(Pcov)';
        
        xEst = xEst + (K * (y - yEst))';
    end
else
    disp('no observation update')

end

%% Global quaternion and reset error GRP
% global quaternion
qeTmp = grp2q(4, fGRP, aGRP, xEst(1,1:3));
qOut = qMult(4, 1, qeTmp, qSigma(1,:));

Pout = Pcov;

% error GRP reset
xEstOut = zeros(1,3);
end
