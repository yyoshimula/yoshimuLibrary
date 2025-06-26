function    [xEstOut, qOut, Pout] = lcCKF(jd, dt_, xEst, Pcov, y, ...
     Qest, Rest, grp, qGlobal, w, r, obsECI, sat, const, earthVSOP)

n_ = length(xEst); % # of state variables

%% initial sigma points
xx = ckfSigma(Pcov, xEst); % 2n x n

% error quaternion sigma points
qeSigma = grp2q(4, grp.f, grp.a, xx(:,1:3)); % 2n x 4

% global quaternion sigma points
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % 2nx4

%[text] ### propagation step
% X0
tmp = qPropMat(4, dt_, w) * qGlobal';  % discrete closed form
qGlobal = tmp';

% X
for j = 1:2*n_
    tmp = qPropMat(4,dt_, w) * qSigma(j,:)';
    qSigma(j,:) = tmp';

    % quaternion to error GRP
    qeSigma(j,:) = qErr(4, qSigma(j,:), qGlobal);
    xx(j,1:3) = q2grp(4, grp.f, grp.a, qeSigma(j,:)); % error GRP sigma point
end
% mean and covariance
xEst = sum(xx,1) / (2 * n_); % predicted state
Pcov = ckfCov(xEst, xx, 1.*Qest);

% sigma point recalculation
xx = ckfSigma(Pcov, xEst);

%% observation step
sunPos = sun(jd, const, earthVSOP); % km, sun position @inertial frame (J2000.0)
sunPos = sunPos .* 10^3; % m
satPos = r; % sat position @inertial frame, m
obsPos = obsECI; % observer position @inertial frame, m

options.BRDF = "AS";
options.mex = "off";

if (y ~= Inf) % observation updateする
    nu = 1; % earth shadowing, とりあえず

    for j = 1:2*n_
        [mApp, ~] = lc(sat, 4, qSigma(j,:), satPos, obsPos, sunPos, nu, options);
        yy(j,1) = mApp;
    end
   
    % predicted mean output 
    yEst = sum(yy, 1) / (2 * n_);

    if(isnan(yEst) || isinf(yEst))
        disp('no observation update');

    else %observation update
        % Calculate correlation
        [Pyy, ~, K] = ckfCorrGain(xEst, xx, yEst, yy, Rest);

        % Update
        Pcov = Pcov - K * Pyy * K';
        
        xEst = xEst + (K * (y - yEst))';
    end
else

    disp('no observation update')

end

%% Global quaternion and reset error GRP
% global quaternion
qeTmp = grp2q(4, grp.f, grp.a, xEst(1,1:3));
qOut = qMult(4, 1, qeTmp, qGlobal);

Pout = Pcov;

% error GRP reset
xEstOut = zeros(1,3);
end


%[appendix]{"version":"1.0"}
%---
