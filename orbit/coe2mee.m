%[text] # calculate modified equinotical elements from classical orbital elements
%[text] modified equinotical elelments: 
%[text] ${\\oe}\_{\\rm MEE} = \[ p, f, g , h, k , L\]  = \\\\\[a(1-e^2),  e\\cos{(w+\\Omega)},  e\\sin{(w+\\Omega)} , \\tan{(i/2)}\\cos{(\\Omega)}, \\tan{(i/2)}\\sin{(\\Omega)}, w+\\Omega+\\nu \]$
%[text] classical orbital elements: 
%[text] ${\\oe} = \[ a, e, i, \\Omega, w , \\nu \~({\\rm or}\~M)\]$
function oe = coe2mee(oe)

oe.p_ = oe.a .* (1 - oe.e.^2);
oe.f_ = oe.e .* cos(oe.w + oe.raan);
oe.g_ = oe.e .* sin(oe.w + oe.raan);
oe.h_ = tan(oe.inc./2) .* cos(oe.raan);
oe.k_ = tan(oe.inc./2) .* sin(oe.raan);
oe.L_ = oe.w + oe.nu + oe.raan;

oe.L_ = mod(oe.L_, 2*pi);

end

%[appendix]{"version":"1.0"}
%---
