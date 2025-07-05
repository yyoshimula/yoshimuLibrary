function ukf_ = setUKFpara(n_, ukf_)
% arguments
%     n_
%     ukf_
% end

%% tuning parameters
if isfield(ukf_, 'alp')
    % as is
else
    ukf_.alp = 1e-4;
end

if isfield(ukf_, 'bet')
    % as is
else
    ukf_.bet = 2;
end

if isfield(ukf_, 'kappa')
    % as is
else
    ukf_.kappa = 3 - n_;
end

if isfield(ukf_, 'lambda')
    % as is
else
    ukf_.lambda = ukf_.alp^2 * (n_ + ukf_.kappa) - n_;
end

%% weights
if isfield(ukf_, 'wm')
    % as is
else
    ukf_.wm(1) = ukf_.lambda / (n_ + ukf_.lambda); % for mean
    ukf_.wm(2:(2*n_+1)) = 1 / (2 * (n_ + ukf_.lambda));
end
ukf_.wm = ukf_.wm(:);

if isfield(ukf_, 'wc')
    % as is
else
    ukf_.wc(1) = ukf_.lambda / (n_ + ukf_.lambda) + 1 - ukf_.alp^2 + ukf_.bet;
    ukf_.wc(2:(2*n_+1)) = ukf_.wm(2:end);
end
ukf_.wc = ukf_.wc(:);

% disp('checking normalization:')
% disp('weights for mean')
% disp(ukf_.w0m + 2 * n_ * ukf_.wim)
% 
% disp('weights for covariance')
% disp(ukf_.w0c + 2 * n_ * ukf_.wic)


end

