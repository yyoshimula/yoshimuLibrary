function    [xEstOut, Pout] = lcUKFonLieGroup(jd, dt_, xEst, Pcov, y, ...
    ukfPara, w, r, obsECI, sat, const, earthVSOP)

% not yet completed

global nWindow adapNoise ind Rest Qest

Qwindow = zeros(1, 1); % adaptive process noise

n_ = length(xEst); % # of state variables
nP_ = length(Pcov);

%% initial sigma points
[x0, xx] = ukfSigmaQ(4, n_, ukfPara.lambda, Pcov, xEst);

%% propagation step
% X0
tmp = qPropMat(4, dt_, w) * x0';  % discrete closed form
x0 = tmp';

% X
for j = 1:2*nP_
    tmp = qPropMat(4, dt_, w) * xx(j,:)';
    xx(j,:) = tmp';
end

% mean and covariance
xEst = x0; % predicted state
% xEst = qAve([x0;xx]);
Pcov = ukfCovQ(4, xEst, x0, xx, ukfPara.w0c, ukfPara.wic, 0 .* Qest);

% sigma point recalculation
[x0, xx] = ukfSigmaQ(4, n_, ukfPara.lambda, Pcov, xEst);

%% observation step
sunPos = sun(jd, const, earthVSOP); % km, sun position @inertial frame (J2000.0)
sunPos = sunPos .* 10^3; % m
satPos = r; % sat position @inertial frame, m
obsPos = obsECI; % observer position @inertial frame, m

yez = zeros(2*nP_, 1);
if (y ~= Inf) % observation updateする
    % nu = 1; % earth shadowing, とりあえず
    nu = shadow(satPos, sunPos, const.RSm, const.REm);    

    [mApp, ~] = lc(sat, 4, x0, satPos, obsPos, sunPos, nu, 'AS');
    ye0 = mApp;

    for j = 1:2*nP_
        [mApp,~] = lc(sat, 4, xx(j,:), satPos, obsPos, sunPos, nu, 'AS');
        yez(j,1) = mApp;
    end    

    % predicted mean output
    ye = ukfPara.w0m .* ye0 + ukfPara.wim .* sum(yez, 1);

    if(isnan(ye) || isinf(ye)) % 予測観測量がinf or NaNになったら観測更新しない
        disp('no observation update due to predicted obs. value');

    else %observation update
        % Calculate correlation and gain
        [Pyy, ~, K] = ukfCorrGainQ(4, xEst, x0, xx, ye, ye0, yez, ukfPara.w0c, ukfPara.wic, Rest);

        % adaptive noise
        dy = y - ye; % output residual
        Qwindow(:,:,ind) = dy' * dy;

        if adapNoise
            sigQ = sum(Qwindow, 3);
            sigQ = sigQ ./ nWindow;
            Qest = K * sigQ * K';
            Qest = nearestSPD(Qest);

            Rest = sigQ;
            % Rest = nearestSPD(Rest);

        else
            % do nothing
        end
        ind = (ind < nWindow) * (ind + 1) + (ind >= nWindow) * 1;

        % Update
        Pcov = Pcov - K * Pyy * K';
        xEst = qMult(4, 1, expQ(4, (K * (y - ye))'), xEst);
        % tmp = expQ(4, (K * (y - ye))');
        % rad2deg(q2zyx(4, tmp))

    end
else
    disp('no observation update')

end

%% Global quaternion and reset error GRP
Pout = Pcov;
xEstOut = xEst;

end
