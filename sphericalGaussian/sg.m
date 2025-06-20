%[text] # Spherical Gaussian
%[text] `x:` x-position on a sphere, NxM matrix
%[text] `y:` y-position on a sphere, NxM matrix
%[text] `z:` z-position on a sphere, NxM matrix
%[text] `p`: lobe directional vector, 1x3 vector 
%[text] `lam`: sharpness, scalar 
%[text] `mu`: amplitude, scalar 
%[text] `G`: spherical Gaussian, Nx1 vector 
%[text] Spherical Gaussian is defined as
%[text] $\\mathcal{G}({\\bf x}, {\\bf p}, \\lambda, \\mu) = \\mu \\exp{(\\lambda({\\bf x}^T{\\bf p} - 1))}$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20200911  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sgMult, sg3.
function G = sg(x, y, z, p, lam, mu)
p = p ./ norm(p);

G = mu .* exp(lam .* (x.*p(1)+y.*p(2)+z.*p(3) - 1));

end

%[appendix]{"version":"1.0"}
%---
