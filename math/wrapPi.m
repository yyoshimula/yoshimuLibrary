%[text] # wrapPi
%[text] 角度 (rad) を半開区間 [-pi, pi) に正規化する（+pi は -pi に写る）
%[text] ## inputs
%[text] - lambda: 角度 [rad]（スカラー・ベクトル・行列可）
%[text] ## outputs
%[text] - lambda: [-pi, pi) に正規化した角度 [rad]
%[text] See also mod
function lambda = wrapPi(lambda)

lambda = mod(lambda + pi, 2*pi) - pi;
end
