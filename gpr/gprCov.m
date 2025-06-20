%[text] # covariance matrix of Gaussian Process Regression
%[text] ## inputs
%[text] `xAst`: $\\bf x^\\ast\n$ test data, m x n matrix 
%[text] `xTrain`: traning data, dxn matrix 
%[text] `yTrain`: traning output data, dx1 vector 
%[text] `Kinv`: inverse matrix of the kernel matrix, dxd matrix 
%[text] `hypPara`: hyper parameters for kernel
%[text] ## outputs
%[text] `yPred`: predicted mean value for test data $\\bf x^\\ast\n$, mxn matrix
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gpCov.
function sigPred = gprCov(xAst, xTrain, Kinv, hypPara)

d = size(xTrain,1); % # of data
m = size(xAst, 1); % # of test data
%[text] $k\_{\\ast \\ast}\n\n$
kAstAst = zeros(m,1);

kAst = zeros(d,m); %行方向に学習データに対する値，列方向に各テストデータに対する値

for i = 1:m
    kAstAst(i,1) = kernelGauss(xAst(i,:), xAst(i,:), hypPara); % mx1

    for j = 1:d
        kAst(j,i) = kernelGauss(xAst(i,:), xTrain(j,:), hypPara); % dxm, test dataが複数なので行列になる
    end   

end

sigPred = kAstAst - diag(kAst' * Kinv * kAst); % mx1

end

%[appendix]{"version":"1.0"}
%---
