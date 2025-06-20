function ukfPara = setUKFpara(n_, ukfPara)
arguments
    n_
    ukfPara
end

%% tuning parameters
if isfield(ukfPara, 'alp')
    % as is
else
    ukfPara.alp = 1e-4;
end

if isfield(ukfPara, 'bet')
    % as is
else
    ukfPara.bet = 2;
end

if isfield(ukfPara, 'kappa')
    % as is
else
    ukfPara.kappa = 3 - n_;
end

if isfield(ukfPara, 'lambda')
    % as is
else
    ukfPara.lambda = ukfPara.alp^2 * (n_ + ukfPara.kappa) - n_;
end

%% weights
if isfield(ukfPara, 'w0m')
    % as is
else
    ukfPara.w0m = ukfPara.lambda / (n_ + ukfPara.lambda); % for mean
end

if isfield(ukfPara, 'wim')
    % as is
else
    ukfPara.wim = 1 / (2 * (n_ + ukfPara.lambda));
end

if isfield(ukfPara, 'w0c')
    % as is
else
    ukfPara.w0c = ukfPara.lambda / (n_ + ukfPara.lambda) + (1 - ukfPara.alp^2 + ukfPara.bet); % for covariance
end

if isfield(ukfPara, 'wic')
    % as is
else
    ukfPara.wic = ukfPara.wim;
end


% disp('checking normalization:')
% disp('weights for mean')
% disp(ukfpara.w0m + 2 * n_ * ukfpara.wim)
% 
% disp('weights for covariance')
% disp(ukfpara.w0c + 2 * n_ * ukfpara.wic)


end

