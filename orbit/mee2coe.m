%[text] # calculate clasical orbital elements from modified equinotical elements
%[text] modified equinotical elelments: 
%[text] ${\\oe}\_{\\rm MEE} = \[ p, f, g , h, k , L\]  = \\\\\[a(1-e^2),  e\\cos{(w+\\Omega)},  e\\sin{(w+\\Omega)} , \\tan{(i/2)}\\cos{(\\Omega)}, \\tan{(i/2)}\\sin{(\\Omega)}, w+\\Omega+\\nu \]$
%[text] classical orbital elements: 
%[text] ${\\oe} = \[ a, e, i, \\Omega, w , \\nu \~({\\rm or}\~M)\]$
function oe = mee2coe(oe)

oe.a = oe.p_ ./ (1 - oe.f_.^2 - oe.g_.^2);
oe.e = sqrt(oe.f_.^2 + oe.g_.^2);
oe.inc = atan2(2*sqrt(oe.h_.^2 + oe.k_.^2), 1 - oe.h_.^2 - oe.k_.^2);
oe.raan = atan2(oe.k_, oe.h_);
oe.w = atan2(oe.g_.*oe.h_ - oe.f_.*oe.k_, oe.f_.*oe.h_ + oe.g_.*oe.k_);
oe.nu = oe.L_ - oe.raan - oe.w;

oe.raan = mod(oe.raan, 2*pi);
oe.w = mod(oe.w, 2*pi);
oe.nu = mod(oe.nu, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
