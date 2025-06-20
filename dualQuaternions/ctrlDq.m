%[text] # control dual quaternions for $C^1$sclerp
%[text] scalar: 0: quaternionのスカラ部q0としてq=\[q0 q1 q2 q3\]という定義
%[text]             4: quaternionのスカラ部q4としてq=\[q1 q2 q3 q4\]という定義
%[text] `dtp`: time step, s 
%[text] `dq1, dq2`: dual quaternions intepolated from dq1 to dq2 (for each) 
%[text] `w1, w2`: corresponding angular rates 
%[text] `v1, v2`: corresponding translational velocities 
%[text] `dqa`, dqb: control dual quaternions, 1x8 (for each)
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
