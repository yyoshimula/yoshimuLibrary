%[text] # Gaussian kernel for each input vector
%[text] ## inputs
%[text] `x, xp`: $\\bf x, x'$ to be calculated, 1xn vector
%[text] `hypPar`: hyper parameters
%[text] ## outputs
%[text] `yPred`: predicted mean value for test data $\\bf x^\\ast\n$, Mxn matrix
%[text] ## note
%[text] $k({\\bf x}\_i, {\\bf x}\_j; \\theta\_1, \\theta\_2) := \\theta\_1 \\exp{\\left( -\\frac{\\|{\\bf x\_i -x\_j} \\|^2}{\\theta\_2}\\right)}$
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230614  y.yoshimura, y.yoshimula@gmail.com
%[text] See also gpCov, kernelGaussMat.
function k = kernelGauss(x, xp, hypPara)
% arguments
%     x (1,:) 
%     xp (1,:) 
%     hypPara (1,2) {mustBeNumeric}    
% end

x = x(:); % nx1 vector
xp = xp(:); % nx1 

k = hypPara(1) * exp(-(x - xp)' * (x - xp) ./ hypPara(2));
  
end

%[appendix]{"version":"1.0"}
%---
