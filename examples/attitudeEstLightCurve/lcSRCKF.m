function    [xEstOut, qOut, Pout] = lcSRCKF(jd, dt_, xEst, Pcov, y, ...
     Qest, Rest, grp, qGlobal, w, r, obsECI, sat, const, earthVSOP)

n_ = length(xEst); % # of state variables

srckf.wm = [2 / (n_ + 2)
    repmat(1/(2*(n_+2)), 2*n_, 1)];  % mean重み
srckf.wc = srckf.wm;

%% initial sigma points
xx = srckfSigma(Pcov, xEst); % (2n+1) x n

% error quaternion sigma points
qeSigma = grp2q(4, grp.f, grp.a, xx(:,1:3)); % (2n+1) x 4

% global quaternion sigma points
qSigma = qMult(4, 1, qeSigma, repmat(qGlobal, length(xx),1)); % (2n+1)x4

%[text] ### propagation step

% X
for j = 1:(2*n_+1)
    tmp = qPropMat(4,dt_, w) * qSigma(j,:)';
    qSigma(j,:) = tmp';

    % quaternion to error GRP
    qeSigma(j,:) = qErr(4, qSigma(j,:), qGlobal);
    xx(j,1:3) = q2grp(4, grp.f, grp.a, qeSigma(j,:)); % error GRP sigma point
end

% mean and covariance
xEst = sum(srckf.wm .* xx, 1); % predicted state
Pcov = srckfCov(xEst, xx, srckf.wc, 1.*Qest);

% sigma point recalculation
xx = srckfSigma(Pcov, xEst);

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
        yy(j,:) = mApp;
    end
   
    % predicted mean output 
    yEst = sum(srckf.wm .* yy, 1);

    if(isnan(yEst) || isinf(yEst))
        disp('no observation update');

    else %observation update
        % Calculate correlation
        [Pyy, ~, K] = srckfCorrGain(xEst, xx, yEst, yy, srckf.wc, Rest);

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
