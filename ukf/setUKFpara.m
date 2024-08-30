function ukfPara = setUKFpara(n_, alp, bet, kap, lam)
arguments
    n_
    alp = 1e-4
    bet = 2
    kap = 3 - n_
    lam = alp^2 + (n_ + kap) - n_    
end

ukfPara.alpha = alp;
ukfPara.beta = bet;

ukfPara.kappa = kap;
ukfPara.lambda = lam;

ukfPara.w0m = ukfPara.lambda / (n_ + ukfPara.lambda); % for mean
ukfPara.wim = 1 / (2 * (n_ + ukfPara.lambda));

ukfPara.w0c = ukfPara.lambda / (n_ + ukfPara.lambda) + (1 - ukfPara.alpha^2 + ukfPara.beta); % for covariance
ukfPara.wic = ukfPara.wim;

end

