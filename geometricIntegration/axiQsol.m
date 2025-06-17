function [qOut, wOut] = axiQsol(t_, qIni, wIni, MOI)
%[text] cf.     Andrle, M. S. & Crassidis, J. L. Geometric Integration of Quaternions. *J. Guid., Control, Dyn.* **36**, 1762–1767 (2013).  
t_ = t_(:);

wn = wIni(3) * (MOI(1,1) - MOI(3,3)) / MOI(1,1);


wOut(:,1) = wIni(1) * cos(wn .* t_) + wIni(2) * sin(wn .* t_);
wOut(:,2) = wIni(2) * cos(wn .* t_) - wIni(1) * sin(wn .* t_);
wOut(:,3) = wIni(3) .* ones(length(t_), 1);

Hini = MOI * wIni(:);
h0 = Hini ./ norm(Hini);
h0 = h0';

wi = norm(Hini) / MOI(1,1);

alp = 0.5 .* wn .* t_;
bet = 0.5 .* wi .* t_;

y = [h0(:,1).*cos(alp).*sin(bet) + h0(:,2).*sin(alp).*sin(bet), ...
    h0(:,2).*cos(alp).*sin(bet) - h0(:,1).*sin(alp).*sin(bet), ...
    h0(:,3).*cos(alp).*sin(bet) + sin(alp).*cos(bet), ...
    cos(alp).*cos(bet) - h0(:,3).*sin(alp).*sin(bet)];

qOut = qMult_(4, 1, y, qIni);


end

%[appendix]{"version":"1.0"}
%---
