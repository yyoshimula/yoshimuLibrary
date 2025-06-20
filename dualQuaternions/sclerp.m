%[text] # Screw linear interpolation of dual quaternions (sclerp)
%[text] `t`: normalzied time, i.e., $0\<t\\leq 1$
%[text] `scalar:` specify the definition of the quaternion 
%[text] `scalar == 0`
%[text] ${\\bf q} = \[q\_0,q\_1,q\_2,q\_3\]^T=\[\\cos(\\theta/2), {\\bf e}^T\\sin(\\theta/2)\]^T$
%[text] `scalar == 4`
%[text] ${\\bf q} = \[q\_1,q\_2,q\_3, q\_4\]^T=\[{\\bf e}^T\\sin(\\theta/2), \\cos(\\theta/2)\]^T$
%[text] `dq1`: starting dual quaternion, 1x8 
%[text] `dq2`: ending dual quaternion, 1x8 
%[text] dqt: interpolated dual quaternions, nx8
function dqt = sclerp(t, scalar, dq1, dq2)
% arguments
%     t (:,1) {mustBeNumeric}
%     scalar (1,1) {mustBeMember(scalar, [0, 4])}
%     dq1 (:,8) {mustBeNumeric}
%     dq2 (:,8) {mustBeNumeric}
% end

t = t(:); % column vector

nT = length(t);

%[text] ## rotation
dqtmp = dqMult(scalar, 0, dqInv(scalar, dq1), dq2);

% scalarに基づく分岐を最初に行い、重複計算を避ける
if scalar == 0
    eTheta = 2 * acos(dqtmp(1));
    eAxis = dqtmp(2:4);
    d = dqtmp(5);
    dqVec = dqtmp(6:8);
else % scalar == 4
    eTheta = 2 * acos(dqtmp(4));
    eAxis = dqtmp(1:3);
    d = dqtmp(8);
    dqVec = dqtmp(5:7);
end

% 共通計算を事前に行う
sinHalfTheta = sin(eTheta/2);
cosHalfTheta = cos(eTheta/2);

% 0除算チェック
if abs(sinHalfTheta) < eps
    % 回転がない場合の処理
    dqt = repmat(dq1, nT, 1);
    return;
end

eAxis = eAxis / sinHalfTheta;

d = -2 * d / sinHalfTheta;
p = (dqVec - d/2 * cosHalfTheta * eAxis) / sinHalfTheta;

% ベクトル化された三角関数計算
tThetaHalf = t * eTheta/2;
cosTthetaHalf = cos(tThetaHalf);
sinTthetaHalf = sin(tThetaHalf);

% 事前に繰り返される値を計算
tdHalf = t * d/2;
eAxisRep = repmat(eAxis', nT, 1);
pRep = repmat(p', nT, 1);

% 三角関数項を事前計算
sinTerm = sinTthetaHalf;
cosTerm = cosTthetaHalf;
neg_t_d_half_sin = -tdHalf .* sinTerm;
t_d_half_cos_eAxis = tdHalf .* cosTerm .* eAxisRep;

% 最終的な dual quaternion の構築
if scalar == 0
    dqInterp = [cosTerm, ...
                eAxisRep .* sinTerm, ...
                neg_t_d_half_sin, ...
                pRep .* sinTerm + t_d_half_cos_eAxis];
else % scalar == 4
    dqInterp = [eAxisRep .* sinTerm, ...
                cosTerm, ...
                pRep .* sinTerm + t_d_half_cos_eAxis, ...
                neg_t_d_half_sin];
end

% 最終的な乗算
dq1rep = repmat(dq1, nT, 1);
dqt = dqMult(scalar, 0, dq1rep, dqInterp);

end

%[appendix]{"version":"1.0"}
%---
