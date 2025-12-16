%[text] # calculate relative position and attitude between chief and deputy
function [chief, deputy, rel] = calcRelPosVelAtti(chief, deputy, anomalyFlag, const)

%[text] ## absolute position and velocity at inertial frame
[chief.rI, chief.vI] = oe2rv(chief.oe, anomalyFlag, const.GE); % km, km/s
[deputy.rI, deputy.vI] = oe2rv(deputy.oe, anomalyFlag, const.GE); % km, km/s, position and vel. at inertial frame
%[text] ### attitude w.r.t. RTN frame (o-frame)
hVec = cross(chief.rI, chief.vI); % orbital angular momentum expressed in inertial frame, nx3
Roi = triad(normRow(chief.rI), normRow(hVec), repmat([1, 0, 0], length(hVec), 1), repmat([0, 0, 1],length(hVec), 1)); % DCM from i-frame to RTN frame, 3x3xn matrix

if isfield(chief, 'q') % if attitude is given
    % DCM
    chief.qoi = zeros(length(hVec), 4); % quaternion from inertial frame to RTN frame
    chief.qbo = zeros(length(hVec), 4); % quaternion from RTN frame to chief's body-fixed frame
    deputy.qbo = zeros(length(hVec), 4); % quaternion from RTN frame to deputy's body-fixed frame
    for i = 1:length(hVec)
        chief.qoi(i,:) = dcm2q(4, Roi(:,:,i));
        chief.qbo(i,:) = qMult(4, 1, chief.q(i,:), qInv(4, chief.qoi(i,:)));
        deputy.qbo(i,:) = qMult(4, 1, deputy.q(i,:), qInv(4, chief.qoi(i,:)));
    end
    % angular rate
    deputy.woi = qRotation(4, repmat([0 0 chief.n],length(hVec),1), deputy.qbo); % i系に対するo系の角速度をdeputyのb系で表したもの
    deputy.wbo = deputy.w - deputy.woi; % deputy's angular rate w.r.t. orbital frame
end

%[text] ## relative position and velocity
% nonlinear relative motion
rel.rNonlinI = deputy.rI - chief.rI; % km, inertial frame
rel.rNonlinRTN = qRotation(4, rel.rNonlinI, chief.qoi); % km, RTN frame
rel.vNonlinI = deputy.vI - chief.vI; % km/s, inertial frame
rel.vNonlinRTN = qRotation(4, rel.vNonlinI, chief.qoi); % km/s, RTN frame
rel.vNonlinRTN = rel.vNonlinRTN - cross(repmat([0, 0, chief.n], size(rel.rNonlinI,1), 1), rel.rNonlinI);

% ROE
rel.roe = oe2roe(chief.oe, deputy.oe, anomalyFlag);
rel.roe(:,2) = wrapPi(rel.roe(:,2)); % wrap angle ( [-\pi, \pi) に変換 )

% mapping from ROE to RTN frame
rel.rMappedRTN = roe2rtn(rel.roe, chief.oe, anomalyFlag, const.GE); % km, RTN frame

end

%[appendix]{"version":"1.0"}
%---
