function flag = selfShadow_obj(sun, n_j, vert_j, n_i, vert_i)
% ----------------------------------------------------------------------
%
%    20200815  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------
% calculate shadow if j-th face make the shadow on the i-th face
% j番目のfaceがi番目のfaceに影を作るかを判定

sun = sun(:);
vert_i = vert_i(:);

a = n_j(1);
b = n_j(2);
c = n_j(3);

d = n_j * vert_j(1,:)'; % 3x1 vector

cos_j = n_j * sun; % scalar
cos_i = n_i * sun;

flag = 1; % nominal

if(cos_j > 0.0 && cos_i > 0.0) 	% どっちも太陽光を受けていたら計算する
    K = d - [a, b, c] * vert_i;
    K = K / ([a, b, c] * sun); % scalar
    
    %     交点計算
    Q = vert_i + K .* sun; % 3x1
    v3 = Q - vert_i;
    inner = v3' * sun;
    
    %	交点がメッシュ内にあるかどうかの計算
    if (inner > 0.0)
        v1 = kron(ones(3,1),Q') - vert_j; % 3x3
        
        %v2 = 3x3 vector 3点を順番につなぐベクトル
        v2 = [vert_j(2,:) - vert_j(1,:)
            vert_j(3,:) - vert_j(2,:)
            vert_j(1,:) - vert_j(3,:)]; % 3x3 matrix
        
        % 外積v2×v1（= crossV）を求める //
        crossV = cross(v2,v1,2); % 3x3
        D = [crossV(1,:) * crossV(2,:)'
            crossV(1,:) * crossV(3,:)'];
        
        if(D(1) >= 0.0 && D(2) >= 0.0)
            flag = 0;          
        else
            flag = 1;
        end
    else
        flag = 1;
    end
else
    %     do nothing    
    flag = 1; % nominal
end

end