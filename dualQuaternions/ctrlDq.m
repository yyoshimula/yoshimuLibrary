%[text] # control dual quaternions for $C^1$ sclerp
%[text] ## inputs
%[text] `dtp`: time step, s
%[text] `dq1`: starting dual quaternion, 1x8
%[text] `dq2`: ending dual quaternion, 1x8
%[text] `w1`: angular rate at dq1, 1x3
%[text] `v1`: translational velocity at dq1, 1x3
%[text] `w2`: angular rate at dq2, 1x3
%[text] `v2`: translational velocity at dq2, 1x3
%[text] ## outputs
%[text] `dqa`: control dual quaternion for the dq1 endpoint, 1x8
%[text] `dqb`: control dual quaternion for the dq2 endpoint, 1x8
%[text] ## note
%[text] NA
%[text] ## references 
%[text] Allmendinger, F., Charaf Eddine, S., & Corves, B. (2018). Coordinate-invariant rigid-body interpolation on a parametric C1 dual quaternion curve. Mechanism and Machine Theory, 121, 731-744. https://doi.org/10.1016/j.mechmachtheory.2017.11.023
%[text] ## revisions
%[text] 20210310  y.yoshimura
%[text] See also sclerp, qMult, dqConj.
function [dqa, dqb] = ctrlDq(dtp, dq1, dq2, w1, v1, w2, v2)

% for dqa
tmp = [dtp/6.*w1, dtp/6.*v1];
eTheta = 2 .* norm(tmp(1:3));
eAxis = tmp(1:3) ./ (eTheta/2);
dqtmp_r = [eAxis.*sin(eTheta/2) cos(eTheta/2)];
dqtmp_d = 0.5 .* qMult(4, 0, [tmp(4:6) 0], dqtmp_r);
dq_tmp = [dqtmp_r dqtmp_d];
dqa = dqMult(4, 0, dq1, dq_tmp);

% for dqb
tmp = [-dtp/6.*w2, -dtp/6.*v2];
eTheta = 2 .* norm(tmp(1:3));
eAxis = tmp(1:3) ./ (eTheta/2);
dqtmp_r = [eAxis.*sin(eTheta/2) cos(eTheta/2)];
dqtmp_d = 0.5 .* qMult(4, 0, [tmp(4:6) 0], dqtmp_r);
dq_tmp = [dqtmp_r dqtmp_d];
dqb = dqMult(4, 0, dq2, dq_tmp);


end

%[appendix]{"version":"1.0"}
%---
