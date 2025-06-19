%[text] # Rotation period of free rotational motion of a rigid body
%[text] ## inputs
%[text] MOI: principal axes of moment of inertia, $\\rm kgm^2$, 3x3 matrix 
%[text] w: angular rate, rad/s, 1x3 vector 
%[text] ## output
%[text] T: period of free rotation, s 
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20190121  y.yoshimura, y.yoshimula@gmail.com
%[text] See also 
function T = rotPeriod(MOI, w)
% arguments
%     MOI (3,3)
%     w (1,3) 
% end

w = w(:); % make column vector

Jx = MOI(1,1);
Jy = MOI(2,2);
Jz = MOI(3,3);

E = 0.5 * w' * MOI * w;
h = MOI * w;
h2 = h' * h;

m = (Jy - Jx) * (2*Jz * E - h2) / (Jz - Jy) / (h2 - 2*Jx*E);
myfun = @(xi) 1 ./ sqrt(1 - m * (sin(xi)).^2);

K = integral(myfun, 0, pi/2);

kappa = sqrt(Jx * (Jz - Jx) / Jy / (Jz - Jy));
w3m = sqrt((h2 - 2*Jx*E) / Jz / (Jz - Jx));

wp = ((Jy - Jz) / Jx) * kappa * w3m;
T = abs(4 * K / wp);

end

%[appendix]{"version":"1.0"}
%---
