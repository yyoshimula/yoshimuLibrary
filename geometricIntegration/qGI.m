%[text] # geometric integration of quaternion (mixed-scheme)
%[text] 
function [qOut, wOut] = qGI(scalar, tspan, qIni, wIni, nGI, MOI)
arguments (Input)
    scalar
    tspan
    qIni
    wIni
    nGI
    MOI
end

arguments (Output)
    qOut
    wOut
end

[a, b, ~] = butcherTable(nGI, 'CG');

if nGI == 3
    nStage = 3;
elseif nGI == 4
    nStage = 5;
end

kw = zeros(nStage, 3);
K = zeros(4,4,nStage);
wi = zeros(nStage,3); %

wi(1,:) = wIni; % intermidiate angular rate

wOut(1,:) = wi(1,:);

% q4 == scalar partとして計算
if scalar == 1
    qOut(1,:) = [qIni(2:4), qIni(1)];
else
    qOut(1,:) = qIni;
end

dt = diff(tspan);

%[text] ## geometric integration
for k = 1:length(tspan)-1

    % when i = 1
    wi(1,:) = wOut(k,:);
    tmp = -MOI^(-1) * cross(wi(1,:)', MOI * wi(1,:)');
    kw(1,:) = tmp';
    K(:,:,1) = 0.5 * qMultMat(4, 1,[wi(1,:), 0]); % 4x4xnStage

    for i = 2:nStage    % for each stage
        tmpw = wi(1,:);
        for j = 1:i-1
            tmpw = tmpw + a(i,j) * dt(k) .* kw(j,:);
        end
        wi(i,:) = tmpw;

        tmp = -MOI^(-1) * cross(wi(i,:)', MOI * wi(i,:)');
        kw(i,:) = tmp';
        K(:,:,i) = 0.5 * qMultMat(4, 1,[wi(i,:), 0]); % 4x4xs
    end

    tmp = zeros(1,3);
    for i =1:nStage
        tmp = tmp + dt(k) .* b(i) .* kw(i,:);
    end
    wOut(k+1,:) = wOut(k,:) + tmp;

    tmp = qOut(k,:)';
    for i = 1:nStage
        % use matrix exponent
        tmp = expm(dt(k) .* b(i) .* K(:,:,i)) * tmp;

        % or use closed form
        % tmp = qExp(4, dt(k)*b(i), wi(i,:)) * tmp;
    end

    qOut(k+1,:) = tmp';

end

end

%[appendix]{"version":"1.0"}
%---
