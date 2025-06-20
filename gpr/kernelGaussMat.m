%[text] # Gaussian kernel matrix for each training data
%[text] ## inputs
%[text] `xTrain: training data`, dxn vector
%[text] `hypPara`: hyper parameters for kernel
%[text] `sigN`: standard deviation of additive noise, scalar
%[text] ## outputs
%[text] `K`: Gaussian kernel matrix, dxd matrix
%[text] ## note
%[text] Statistics and Machine Learning Toolbox is required
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gpCov.
function K = kernelGaussMat(xTrain, hypPara, sigN)
% arguments
%     xTrain 
%     hypPara (1,2) {mustBeNumeric}
%     sigN (1,1) {mustBeNumeric}
% end

d = gather(size(xTrain,1)); % data length

%[text] #### (Statistics and Machine Learning Toolbox is required)
K = hypPara(1) .* exp(-squareform(pdist(xTrain)).^2 ./ hypPara(2));

%[text] #### use below, if you do not have Statistics and Machine Learning Toolbox
% d = gather(size(xTrain,1)); % data length
% K = zeros(d,d);
% for i = 1:d 
%     for j = i:d % 上三角部分だけ計算
%         K(i,j) = kernelGauss(xTrain(i,:), xTrain(j,:), hypPara);
%         K(j,i) = K(i,j);
%     end
% end

%[text] #### add noise term
K = K + sigN^2 .* eye(d);

end

%[appendix]{"version":"1.0"}
%---
