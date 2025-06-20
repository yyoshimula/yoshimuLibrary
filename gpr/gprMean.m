%[text] # mean value of Gaussian Process Regression
%[text] ## inputs
%[text] `xAst`: $\\bf x^\\ast\n$ test data, m x n matrix 
%[text] `xTrain`: traning data, Nxn matrix 
%[text] `yTrain`: traning output data, Nxm vector 
%[text] `Kinv`: inverse matrix of the kernel matrix, NxN matrix 
%[text] `hypPara`: hyper parameters for kernel
%[text] ## outputs
%[text] `yPred`: predicted mean value for test data $\\bf x^\\ast\n$, Mxn matrix
%[text] logP: log marginal likelihood, $\\log{p({\\bf y}|X)}$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gpCov.
function [yPred, logP] = gprMean(xAst, xTrain, xMean, yTrain, yMean, L, hypPara)
% arguments
%     xAst
%     xTrain
%     xMean
%     yTrain
%     yMean
%     L
%     hypPara
% end

%[text] ### dimension
d = size(xTrain, 1); % # of data
m = size(xAst, 1); % # of test data

%[text] ## calculating $\\bf k\_\\ast$
kAst = zeros(d,m); %行方向に学習データに対する値，列方向に各テストデータに対する値

for i = 1:m
    for j = 1:d
        kAst(j,i) = kernelGauss(xAst(i,:), xTrain(j,:), hypPara); % dxm matrix
    end
end
%[text] ## mean
%[text] #### directly inverse the kernel matrix
% yPred = xMean + kAst' * Kinv * (yTrain - xMean); % Mxm matrix, m is # of output variables

%[text] #### use the Cholesky decomposition (numerically stable and fast)
alp = L' \ (L \ (yTrain - xMean));
 
yPred = xMean + kAst' * alp; % Mxm matrix, m is # of output variables

%[text] ## log marginal likelihood
logP = -0.5 * yTrain' * alp - sum(diag(L)) - d / 2 * log(2 * pi);

end

%[appendix]{"version":"1.0"}
%---
