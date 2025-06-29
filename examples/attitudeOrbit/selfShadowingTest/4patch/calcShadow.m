% 全ての機器，メッシュに対して処理を行う
% //■影の処理関数(nun.1が影を作るほう，nun.2が影が写るほうのメッシュ)
function face = calcShadow(face, vert, sun)
global MAX_F MAX_EQUIP EQUIP_NUM

% 	// 平面の方程式を求める
% 	// 平面上の点P1(px1,py1,pz1), 点P2(px2,py2,pz2), 点P3(px3,py3,pz3)
% 	// 平面の方程式 a*x + b*y + c*z = d; ・・・(1)
% 	// また、P1P2 = u = (u1,u2,u3), P1P3 = v = (v1,v2,v3) とすると、
% 	// 平面の方程式で下記の式が成立。
% 	// |u2 v2| |u3 v3| |u1 v1|
% 	// |u3 v3| * (x-x0) + |u1 v1|*(y-y0) + |u2 v2|*(z-z0) = 0 ・・・（2）


% 	double Inner;			// 内積計算
% 	double Inner3;
% 	double Q[3];			// 交点
% 	double v1[4][3];		// 4点→交点のベクトル
% 	double v2[4][3];		// 0→1,1→2,2→3,3→0のベクトル
% 	double cross[4][3];	// V2×V1
% 	double K;				// 交点計算用係数その1

% 	//半直線判別関係
% 	double v3[3];			// 中心点→交点ベクトル
% 	double Inner2;		// 内積：v3*sun_light

sun = sun(:); % column vector

%% 係数処理
for m = 1:MAX_EQUIP
    for k = 1:EQUIP_NUM(m)
        n_ = length(face(m,k).pos);
        % 	// 平面方程式の係数 全てのメッシュに対して係数を算出しておく
        face(m,k).a(:,1) = face(m,k).normal(1,:)'; % nx1 vector
        face(m,k).b(:,1) = face(m,k).normal(2,:)';
        face(m,k).c(:,1) = face(m,k).normal(3,:)';
        
        face(m,k).d(:,1) = face(m,k).a .* vert(m,k).coord(face(m,k).make(:,1),1)...
            + face(m,k).b .* vert(m,k).coord(face(m,k).make(:,1),2)...
            + face(m,k).c .* vert(m,k).coord(face(m,k).make(:,1),3); %仮に面を作る4点の1番目の点を使う．nx1 vector
        
        %j番目の太陽方向ベクトルに対して計算を行う
        face(m,k).inner(:,1) = face(m,k).normal' * sun; %全メッシュの内積　nx1 vector
        
        face(m,k).flag = ones(n_,1); % 影判定用フラグ 0のとき，影になっているとする
    end
end

%% self-shadowing
for m = 1:MAX_EQUIP % each equipment
    for k = 1:EQUIP_NUM(m) % each face of an equipment
        for i = 1:MAX_F(m,k) % i番目のfacetが影になるかを計算, each facet
            
            temp = [vert(m,k).coord(face(m,k).make(i,1),:)
                vert(m,k).coord(face(m,k).make(i,2),:)
                vert(m,k).coord(face(m,k).make(i,3),:)
                vert(m,k).coord(face(m,k).make(i,4),:)]; % 4x3
            center = sum(temp,1) ./ 4.0; % the center of i-th facet, 1x3 vector
            
            % 自身の機器からの影
            for kk = 1:EQUIP_NUM(m) % each face
                for ll = 1:MAX_F(m,kk) % each facet
                    coord_ll = [vert(m,kk).coord(face(m,kk).make(ll,1),:)
                        vert(m,kk).coord(face(m,kk).make(ll,2),:)
                        vert(m,kk).coord(face(m,kk).make(ll,3),:)
                        vert(m,kk).coord(face(m,kk).make(ll,4),:)]; % 4x3
                    
                    flag = selfShadow(sun, face(m,kk).normal(:,ll), coord_ll, face(m,k).normal(:,i), center);
                    face(m,k).flag(i) = face(m,k).flag(i) * flag; % OR演算っぽく
                end
            end
            
            % 影の計算，他機器→自身の影 //
            for mm = 1:m-1
                for kk = 1:EQUIP_NUM(mm)
                    for ll = 1:MAX_F(mm,kk)
                        coord_ll = [vert(mm,kk).coord(face(mm,kk).make(ll,1),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,2),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,3),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,4),:)]; % 4x3
                        
                        flag = selfShadow(sun, face(mm,kk).normal(:,ll), coord_ll, face(m,k).normal(:,i), center);
                        face(m,k).flag(i) = face(m,k).flag(i) * flag; % OR演算っぽく
                    end
                end
            end
            
            for mm = m+1:MAX_EQUIP
                for kk = 1:EQUIP_NUM(mm)
                    for ll = 1:MAX_F(mm,kk)
                        coord_ll = [vert(mm,kk).coord(face(mm,kk).make(ll,1),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,2),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,3),:)
                            vert(mm,kk).coord(face(mm,kk).make(ll,4),:)]; % 4x3
                        
                        flag = selfShadow(sun, face(mm,kk).normal(:,ll), coord_ll, face(m,k).normal(:,i), center);
                        face(m,k).flag(i) = face(m,k).flag(i) * flag; % OR演算っぽく

                    end
                end
            end
            
        end
    end
end
end