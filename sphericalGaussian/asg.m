%[text] # Anisotropic Spherical Gaussian
%[text] `vx, vy, vz:` position on sphere, NxM matrix
%[text] `x`, tangent axis, 1x3 vector
%[text] `y`, bi-tangent axis, 1x3 vector
%[text] `z`, lobe axis, 1x3 vector
%[text] `lam`, bandwidth along x-axis, scalar
%[text] `mu`, bandwidth along y-axis, scalar
%[text] `c`, amplitude, scalar
%[text] `asg`, anisotropic spherical Gaussian, Nx1 vector
%[text] 
%[text] Anisotropic spherical Gaussian is defined as
%[text] $\\mathcal{G}({\\bf x}, {\\bf p}, \\lambda, \\mu) = \\mu \\exp{(\\lambda({\\bf x}^T{\\bf p} - 1))}$
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20231004  y.yoshimura, y.yoshimula@gmail.com
%[text] See also sgMult, sg.
function aSG = asg(vx, vy, vz, x, y, z, lam, mu, c)


VX = vx .* x(1) + vy .* x(2) + vz .* x(3);
VY = vx .* y(1) + vy .* y(2) + vz .* y(3);
VZ = vx .* z(1) + vy .* z(2) + vz .* z(3);

S = max(0, VZ); % Nx1 vector
aSG = c .* S .* exp(-lam .* VX.^2 - mu .* VY.^2);


end

%[appendix]{"version":"1.0"}
%---
