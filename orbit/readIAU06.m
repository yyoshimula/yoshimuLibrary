%[text] # read IAU06 values for calculating precession and nutation based on IAU-2006/2000 theory
%[text] ## note
%[text] `iau06x.dat, iau06y.dat, iau06z.dat`
%[text] are required.
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.213
%[text] ## revisions
%[text] 20230608  y.yoshimura, y.yoshimula@gmail.com
%[text] See also leapS
function iau06 = readIAU06(~)

%[text] ## read coefficients
% xys values
iau06x = readmatrix('iau06x.dat');
iau06.axs0 = iau06x(:,2:3);  % reals
iau06.a0xi = iau06x(:,4:17); % integers

iau06.axs0= 10^(-6) * arcs2rad(iau06.axs0);  % micro arcs to rad

iau06y = readmatrix('iau06y.dat');
iau06.ays0 = iau06y(:,2:3);
iau06.a0yi = iau06y(:,4:17);

iau06.ays0 = 10^(-6) .* arcs2rad(iau06.ays0);  % micro arcs to rad

iau06s = readmatrix('iau06s.dat');
iau06.ass0 = iau06s(:,2:3);
iau06.a0si = iau06s(:,4:17);

iau06.ass0 = 10^(-6) * arcs2rad(iau06.ass0);

end



%[appendix]{"version":"1.0"}
%---
