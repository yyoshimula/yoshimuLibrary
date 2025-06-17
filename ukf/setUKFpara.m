function ukfpara = setUKFpara(n_, ukfpara)
arguments
    n_
    ukfpara
end

%% tuning parameters
if isfield(ukfpara, 'alp')
    % as is
else
    ukfpara.alp = 1e-4;
end

if isfield(ukfpara, 'bet')
    % as is
else
    ukfpara.bet = 2;
end

if isfield(ukfpara, 'kappa')
    % as is
else
    ukfpara.kappa = 3 - n_;
end

if isfield(ukfpara, 'lambda')
    % as is
else
    ukfpara.lambda = ukfpara.alp^2 * (n_ + ukfpara.kappa) - n_;
end

%% weights
if isfield(ukfpara, 'w0m')
    % as is
else
    ukfpara.w0m = ukfpara.lambda / (n_ + ukfpara.lambda); % for mean
end

if isfield(ukfpara, 'wim')
    % as is
else
    ukfpara.wim = 1 / (2 * (n_ + ukfpara.lambda));
end

if isfield(ukfpara, 'w0c')
    % as is
else
    ukfpara.w0c = ukfpara.lambda / (n_ + ukfpara.lambda) + (1 - ukfpara.alp^2 + ukfpara.bet); % for covariance
end

if isfield(ukfpara, 'wic')
    % as is
else
    ukfpara.wic = ukfpara.wim;
end


% disp('checking normalization:')
% disp('weights for mean')
% disp(ukfpara.w0m + 2 * n_ * ukfpara.wim)
% 
% disp('weights for covariance')
% disp(ukfpara.w0c + 2 * n_ * ukfpara.wic)


end

