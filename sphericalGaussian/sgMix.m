%[text] # Spherical Gaussian mixture (SG multiplication)
%[text] spherical Gaussian:
%[text] $\\mathcal{G}({\\bf x}, {\\bf p}, \\lambda, \\mu) = \\mu \\exp{(\\lambda({\\bf x}^T{\\bf p} - 1))}$
%[text] spherical Gaussian mixture:
%[text] $\\mathcal{G}({\\bf{x}};{\\bf{p}}\_{1},\\lambda\_{1}, \\mu\_{1})\\mathcal{G}({\\bf{x}};{\\bf{p}}\_{2},\\lambda\_{2}, \\mu\_{2}) = \\mathcal{G}({\\bf{x}}; \\bar{{\\bf{p}}}\_{m}, \\lambda\_{m}, \\mu\_{m})$
%[text] where
%[text] ${\\bf{p}}\_{m} = \\lambda\_{1}{\\bf{p}}\_{1}+\\lambda\_{2}{\\bf{p}}\_{2} \\\\\n\\lambda\_{m} = \\|{\\bf{p}}\_{m}\\| \\\\\n\\bar{{\\bf{p}}}\_{m} = \\frac{{\\bf{p}}\_{m}}{\\|{\\bf{p}}\_{m}\\|} \\\\\n\\mu\_{m} = \\mu\_{1}\\mu\_{2}e^{(\\lambda\_{m}-(\\lambda\_{1} + \\lambda\_{2}))}$
%[text] ## inputs & outputs
%[text] `p1, p2, p3:` lobe axes, nx3 vector
%[text] `lam1, lam2, lam3:` sharpness, scalar,
%[text] `mu1, mu2, mu3:` amplitudes of the SGs
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20210206  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sg.
function [p3, lam3, mu3] = sgMix(p1, lam1, mu1, p2, lam2, mu2)
% arguments
%     p1 (:,3) {mustBeNumeric}
%     lam1 (:,1) {mustBeNumeric}
%     mu1 (:,1) {mustBeNumeric}
%     p2 (:,3) {mustBeNumeric}
%     lam2 (:,1) {mustBeNumeric}
%     mu2 (:,1) {mustBeNumeric}
% end

%[text] ### sharpness
lam3 = vecnorm(lam1.*p1 + lam2.*p2, 2, 2); % nx1

%[text] ### lobe axis
p3 = (lam1.*p1 + lam2.*p2) ./ lam3;

%[text] ### amplitude
mu3 = mu1 .* mu2 .* exp(lam3 - (lam1 + lam2));

end

%[appendix]{"version":"1.0"}
%---
