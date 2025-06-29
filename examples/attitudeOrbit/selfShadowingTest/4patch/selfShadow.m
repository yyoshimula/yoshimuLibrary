function flag = selfShadow(sun, n_1, coord_1, n_2, coord_2)
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

coord_2 = coord_2(:);

a = n_1(1);
b = n_1(2);
c = n_1(3);

d = n_1' * coord_1(1,:)';

cos1 = sun' * n_1;
cos2 = sun' * n_2;

flag = 1; % nominal

if(cos1 > 0.0 && cos2 > 0.0) 	% どっちも太陽光を受けていたら計算する
    K = d - [a, b, c] * coord_2;
    K = K / ([a, b, c] * sun); % scalar
    
    %     交点計算
    Q = coord_2 + K .* sun; % 3x1
    v3 = Q - coord_2;
    inner = v3' * sun;
    
    %	交点がメッシュ内にあるかどうかの計算
    if (inner > 0.0)
        
        v1 = kron(ones(4,1),Q') - coord_1; % 4x3
        
        %v2 = 4x3 vector 4点を順番につなぐベクトル
        v2 = [coord_1(2,:) - coord_1(1,:)
            coord_1(3,:) - coord_1(2,:)
            coord_1(4,:) - coord_1(3,:)
            coord_1(1,:) - coord_1(4,:)];
        
        % 外積v2×v1（= crossV）を求める //
        crossV = cross(v2,v1);
        
        D = [dot(crossV(1,:), crossV(1,:))
            dot(crossV(1,:), crossV(3,:))
            dot(crossV(1,:), crossV(4,:))]; % 3x1
        
%        内外判定
        %         if (D(1) > 0.0 && D(2) > 0.0 && D(3) > 0.0) % こっちでもok
        if((crossV(1,3) > 0.0 && crossV(2,3) > 0.0 && crossV(3,3) > 0.0 && crossV(4,3) > 0.0 ) || ...
                (crossV(1,3) < 0.0 && crossV(2,3) < 0.0 && crossV(3,3) < 0.0 && crossV(4,3) > 0.0 ))
            flag = 0;
        else
            flag = 1;
        end
    else
        flag = 1;
    end
else
    %     do nothing
    
end


end