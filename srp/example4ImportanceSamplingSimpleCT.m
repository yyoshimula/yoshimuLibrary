%[text] # example for importance sampling of rendering equation (Cook–Torrance model)
%[text] モンテカルロ積分に使う乱数の数は `nMCint`
%[text] `それを nMC` 回行って平均と分散を計算し，plotする．（ややこしい）
%[text] 一応，SRPの計算にも対応（ただし太陽定数とかをかけていない）
clc
clear
cls

format long
%[text] ## Monte-Carlo condition & pre-allocation
nMC = 100;

srpTrue = zeros(nMC, 3);
srpImp = zeros(nMC, 3);
%[text] ## Cook–Torrance model
sat.F0 = 0.5 %[output:3a3132f0]
sat.mCT = 0.3;

%[text] ### normal vector
% sat.normal = [ 0 0 1 ]; 
sat.normal = rand(1,3);
sat.normal = sat.normal ./ vecnorm(sat.normal,2,2) %[output:897f66e1]

sat.rho = 0.5;
sat.area = 1;

thetaI = deg2rad(60);
phiI = deg2rad(80);

sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)];
sunB = sunB ./ norm(sunB);

% sunB = [ -0.573721081964441   0.434806422007045   0.694109138025847]

% Monte-Carlo integration condition
nMCintUni = [100, 500, 1000, 5000, 10000, 50000, 100000];

tic
for i = 1:length(nMCintUni)
    N = nMCintUni(i);
    for j = 1:nMC
%[text] ## speuclar only, uniform distribution
        % sampling for sphere region(全球) (note: not hemisphere!)
        thetaR = pi .* rand(N,1); % Nx1
        phiR = 2 * pi .* rand(N,1);

        v = [sin(thetaR).*cos(phiR) sin(thetaR).*sin(phiR) cos(thetaR)]; %Nx3       
        h = sunB + v;
        h = h ./ vecnorm(h,2,2); % bisector vector, Nx3
        thetaH = acos(h(:,3)); % Nx1
        thetaH = acos(sat.normal*h')';        

        VH = dot(v, h, 2)'; % 1xN
        NV = sat.normal * v'; % nxN, 列が各reference vector, vの値
        NS = sat.normal * sunB'; % nx1
        NH = sat.normal * h'; % nxN
        D = exp(-(tan(thetaH')./sat.mCT).^2); % Beckmann distribution, nxN        
        D = D ./ pi ./ sat.mCT.^2 ./ cos(thetaH').^4;        
        %     D = exp(-(thetaH ./ sat.mCT).^2); % Gaussian distribution, nxN

        nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0)); % nx1
        g = sqrt(nest.^2 + VH.^2 - 1); % nxN

        temp1 = 2 * NH .* NV ./ VH; % nxN
        temp2 = 2 * NH .* NS ./ VH; % nxN
        G = min(1, temp1);
        G = min(G, temp2); % nxN

        temp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
        temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
        F = temp1 .* temp2; % nxN

        temp = D .* G .* F ./ NS ./ NV / 4; % nxN
        temp = temp .* NV .* sin(thetaR)'; % nxN

        % for SRP
        csCTx = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,1)'; % nxN matrix
        csCTy = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,2)'; % nxN matrix
        csCTz = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,3)'; % nxN matrix

        % % for rendering eq.
        csCTx = (NS > 0) .* (NV > 0) .* temp;
        csCTy = (NS > 0) .* (NV > 0) .* temp;
        csCTz = (NS > 0) .* (NV > 0) .* temp;

        % probability
        pCT = 1 / pi * 1 / (2 * pi); % scalar

        % specular component, nx3
        srpCs = [sum(csCTx, 2), sum(csCTy, 2), sum(csCTz, 2)];
        srpCs = srpCs ./ pCT ./ N;

        srpTrue(j,:) = srpCs;
    end
    stdTrue(i,:) = std(srpTrue);
    meanTrue(i,:) = mean(srpTrue);
end
toc %[output:2a5a643f]
%[text] ## 
%[text] ## importance sampling
% Monte-Carlo integration condition
nMCintImp = [100 10^4];
nMCintImp = nMCintUni;

eAngle = acos(sat.normal(:,3));
eAxis = cross([0 0 1], sat.normal) .* (eAngle > 1e-5) ...
     + [0 0 1] * (eAngle < 1e-5);
eAxis = eAxis ./ vecnorm(eAxis, 2, 2);

qlb = [eAxis.*sin(eAngle/2), cos(eAngle/2)];

sLocal = qRotation(4, sunB, qlb);
nLocal = qRotation(4, sat.normal, qlb) %[output:9b739829]


tic
for i = 1:length(nMCintImp)
    N = nMCintImp(i);
    for j = 1:nMC
        u1 = rand(N,1);

        thetaH = atan(sqrt(-sat.mCT.^2 .* log(u1))); % Nx1
        phiH = 2 * pi .* rand(N,1);

        h = [sin(thetaH).*cos(phiH) sin(thetaH).*sin(phiH) cos(thetaH)]; %Nx3
        % h = qRotation(4, h, repmat(qlb,N,1));
        v = 2 * (h * sLocal') .* h - sLocal; % Nx3

        VH = dot(v, h, 2)'; % 1xN
        NV = nLocal * v'; % nxN, 列方向にreference vector, vの値
        NS = nLocal * sLocal'; % nx1
        NH = nLocal * h'; % nxN
        SH = (h * sLocal')'; % 1xN

        nest = (1 + sqrt(sat.F0)) ./ (1 - sqrt(sat.F0)); % nx1
        g = sqrt(nest.^2 + VH.^2 - 1); % nxN

        temp1 = 2 * NH .* NV ./ VH; % nxN
        temp2 = 2 * NH .* NS ./ VH; % nxN
        G = min(1, temp1);
        G = min(G, temp2); % nxN

        temp1 = (g - VH).^2 / 2 ./ (g + VH).^2;
        temp2 = (1 + (VH .* (g + VH) - 1).^2 ./ (VH .* (g - VH) + 1).^2);
        F = temp1 .* temp2; % nxN

        % weight
        W = abs(SH) .* G ./ NS ./ NH;

        temp = W .* F ; % nxN

        % for SRP
        % csCTx = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,1)'; % nxN matrix
        % csCTy = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,2)'; % nxN matrix
        % csCTz = NS .* (NS > 0) .* (NV > 0) .* temp .* v(:,3)'; % nxN matrix
        % srpCs = [mean(csCTx), mean(csCTy), mean(csCTz)];
        % srpCs = qRotation(4, srpCs, qInv(4, qlb));        

        % for rendering eq.
        csCTx = (NS > 0) .* (NV > 0) .* temp ;
        csCTy = (NS > 0) .* (NV > 0) .* temp ;
        csCTz = (NS > 0) .* (NV > 0) .* temp ;
        srpCs = [mean(csCTx), mean(csCTy), mean(csCTz)];
        

        srpImp(j,:) = srpCs;
    end
    stdImp(i,:) = std(srpImp);
    meanImp(i,:) = mean(srpImp);
end

toc %[output:2ed18787]

%[text] ## show FIgs
figure %[output:1cedac6f]
tiledlayout(3,1), nexttile %[output:1cedac6f]
plotStd(nMCintUni, meanTrue(:,1), stdTrue(:,1), 'r'), hold on %[output:1cedac6f]
plotStd(nMCintImp, meanImp(:,1), stdImp(:,1), 'b', 'b'); %[output:1cedac6f]


nexttile %[output:1cedac6f]
plotStd(nMCintUni, meanTrue(:,2), stdTrue(:,2)); hold on %[output:1cedac6f]
plotStd(nMCintImp, meanImp(:,2), stdImp(:,2),'b', 'b'); %[output:1cedac6f]

nexttile %[output:1cedac6f]
plotStd(nMCintUni, meanTrue(:,3), stdTrue(:,3)); hold on %[output:1cedac6f]
plotStd(nMCintImp, meanImp(:,3), stdImp(:,3), 'b', 'b'); %[output:1cedac6f]

%[text] ### mean and standard deviation
[meanTrue(end,:) %[output:group:1775564d] %[output:14f0391e]
    meanImp(end,:)] %[output:group:1775564d] %[output:14f0391e]

[norm(meanTrue(end,:)) %[output:group:48f280a3] %[output:4aa62d23]
    norm(meanImp(end,:))] %[output:group:48f280a3] %[output:4aa62d23]
%[text] ## 
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40.7}
%---
%[output:3a3132f0]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"sat","value":"    F0: 0.500000000000000\n"}}
%---
%[output:897f66e1]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"sat","value":"        F0: 0.500000000000000\n       mCT: 0.300000000000000\n    normal: [0.696409420356347 0.512186699201015 0.502671567129578]\n"}}
%---
%[output:2a5a643f]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 1.889894 秒です。\n","truncated":false}}
%---
%[output:9b739829]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"nLocal","rows":1,"type":"double","value":[["0.000000000000000","0","1.000000000000000"]]}}
%---
%[output:2ed18787]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 1.152523 秒です。\n","truncated":false}}
%---
%[output:1cedac6f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAmgAAAHOCAYAAAAsQvUnAAAAAXNSR0IArs4c6QAAIABJREFUeF7t3QuwVWX9\/\/Evt4i4RWiU4hFFKdEYHEpASBFi0F8F\/JEcZOSPZeEFjg6QCJmRGiAlAXEJJAe5GYEoJaBABiReYCgSEUOCMX8dUcNhkCEGCPnP56m9\/\/sc9jln7732Pns9+3mvmaY6rLX2s17ru9f67Gc9a616Z86cOWNMCCCAAAIIIIAAArERqEdAi82+oCEIIIAAAggggIATIKBRCAgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBRAwgggAACCCCAQMwECGgx2yE0BwEEEEAAAQQQIKBFrIFhw4bZ9u3bI66FxRFAAAEEEECgNoGrrrrKlixZUttsJfHvBLSIu\/Gcc86x8vJymzhxYsQ1FWfxL3zhC7Z3797ifHjET\/W57dp02h+xACIsjn0EvDwsin8eEHNcBfY5whVhsWAC2smTJ23lypW2e\/du69q1q\/Xr18+aNGlyFrnmW7Fihb322mt2\/vnn280332xt2rSpdtcQ0IpQtf\/9SA40xbP3PWBSO9ROFAGf68fntvt+3Mm25oIJaKNHj7YjR47Y8OHDbdeuXS6Ebdiw4ayQdvvtt1uLFi1s0KBBtmbNGnv11Vft+eeft0aNGqW1LSsrs\/bt29umTZuytY\/F\/D5\/WX1ueykcaHz297nt1E7xD50+14\/PbS+F2s+meoMIaAcOHLChQ4fali1brHHjxs5n5MiR1rdvXxs4cGDS64033rC7777b1q1b5+Y7c+aMLVu2zM3TrFmztK4dOnSwffv2mf473TRq1Ch3CTSuk89fVp\/bXgoHGp\/9fW47tVP8o6nP9RP3ts+ePdtmzZpV4072dVhOtpUbREBbtWqVbd261aZPn570Wb16te3fv9\/Gjh2b\/NuCBQvs9OnTpoH\/O3bssJYtW1rnzp1rNE0ENIU5Hyd9GRQifZx8bru8aX\/xqg774tlT+9hHEYh7wIyybVWXDSKgzZs3z06dOlWpJ2vbtm22aNEimzt3btLkwQcfdGPUNG+PHj1M81xwwQU2bdo0q1+\/frU9aBqHNnnyZOvVq1c+9w3rQgABBBBAAIEUAQJaiZXD0qVL7fDhw5UC2saNG029aHPmzElu7cMPP+x62jTmrF69eu4S59e+9jXX89apU6dqA1r\/\/v2tefPm3t7JWWK7m81BAAEEEChRAQJaie1Y9YTpuSm6rJGY5s+f724aGDduXPJvCnJvvvmmTZo0Kfm3CRMmuLs+U8eqpfLoEicBrcQKhs1BAAEEEIilAAEtlrsl90YdPXrU3RCgS50aU1ZRUWFDhgwxjTlTwNq5c6d17NjRPvzwQxe2dGPAZZdd5ua78cYb3f\/XnZrpJgJa7vuFJRFAAAEEEMhGgICWjZYn8+rRGgsXLrQ9e\/a4nrPx48e7XrHjx4+78WZPP\/20tWvXznTp86GHHnJbpbFo6kEbMGBAtVupgDZlyhTXO+frozY82YU0EwEEEEAgcAECWgkXgAJZugfUpm6yxp4dO3as2kdrpM6rgLZ582b3UFtf7+Qs4d3NpiGAAAIIlJAAAa2EdmahNyU1oCV60Libs9DqrB8BBBBAIEQBAlqIez3HbU4ENC0+ePBgt5ZXXnnFPXIj8X5OAluOuCyGAAIIIIBAigABjXLIWCA1oKUupJCm56cR2DKmZEYEEEAAAQRqFCCgUSAZC1QX0KquIDWwJV43RQ9bxszMiAACCCCAgBHQKIKMBTINaFV71\/T\/Ez1sBLaMuZkRAQQQQCBgAQJawDs\/203PJaCl610jsGUrz\/wIIIAAAqEJENBC2+MRtjcfAa22wKbLo1Wn6m48qPr3a6+9Nu3WceNChJ3OoggggAACRREgoBWF3c8PLURAy0QiXWjTclX\/\/vLLL6ddXXXL1\/bZdRHsCvkZ1QXW2rY7m38vZPuzaQfzIoAAAqUmQEArtT1awO0pKyuztm3b2tixY6179+4F\/KR4rDrXYJdN6wv5GdUF1mzaV9u8hWy\/PrsuAmChP6PQQbnQ7a+tBvh3BBAojAABrTCuJbnWSy+9wvr0mWrz5z9q3bufcNt4ddVLkpkEt6uvztinpGJgJjYZy4QxY6EDoBQLHWTrYhvCqAa\/t5Ig7ff+K0brX3\/9dTt06FAxPrrOP7PeGd5PFAm9rGyYtWnzuO3YMcvMtvx3XZszWOfJDOapbpb\/BEEmBBBAAAEEQhL41Kc+5V7FGMJEQIu4l\/9P2SXW5n+PWzP70GzsqIzXls\/LoflcV8YbUEczlkJPSylsQ130qtVFSZXKvqgLqxA\/gx69+O91etDiv49i00KNQbvppptsjC7VFeFyXboTTrq\/pbtk5cPJqhQOmKWwDfrCFXrcWF18qUtlX9SFFZ+BQBwFGIMWx70S0zYV4i5OBadEeEoNVpk+biPdSSjdyZWTVUyLimYhgAACCKQVIKBRGBkLRA1oVcNY4kXrifCUGqwIVBnvFmZEAAEEEChBAQJaCe7UQm1SNgEtkzBGCCvUnmK9CCCAAAK+CxDQfN+Dddj+6gIaYawOdwIfhQACCCAQhAABLYjdnJ+NrBrQFMwGDx7sHiaaepmSnrH8eLMWBBBAAIFwBQho4e77rLe8akBTOJs8eXKdPO0968ayAAIIIIAAAh4LENA83nnVNf3999+3Z555xioqKuz666+3bt26WYMGDWrc0sWLF9s111xj7dq1q3a+1ICm3rPZs2fbpk2bSlCQTUIAAQQQQKC4AgS04vrn\/dOPHDninlXWs2dP6927t61YscLq169v06dPr\/azVq9ebffdd589+eST1qVLl4wCGr1ned91rBABBBBAAIGkAAGtxIph5cqVtnXrVps5c6bbshMnTriHbj711FPuRedVp71799rw4cPtoosuciGtc+fOtQY0es9KrGjYHAQQQACB2AkQ0GK3S6I1aPz48S5kDRkyJLmiCRMmuEudVR\/gevToURs4cKCNGzfO1q1bZ7fddpt16tSpxoB2+vRpO3jwoJ1zzjnWpEmTSvOOGjXKysvLo20ASyOAAAIIIBCAgIYJzZqld1tXP6kTJYQpiHdxqjfsjjvusNR3VqoAWrRo4XrKEpPeG69Adfnll9tdd93l\/veIESNqDWibN2+2888\/33jvfAhfGbYRAQQQQKBYAvSgFUu+QJ97zz33uN6z1IA2adIku\/DCC+2WW25JfuqcOXNs\/\/79Nm3aNKtXr14yoHXs2NEaNmyYtnW6SWDKlCncHFCgfcdqEUAAAQQQSAgQ0EqsFhS8WrZsWSmM3Xrrra737Lrrrktubb9+\/dzNA40aNXJ\/Uzdqq1atXA\/aY489RkArsbpgcxBAAAEE\/BIgoPm1v2pt7c6dO03j0HRHZuvWrU0vIB8zZoxt2bLF3TDw9ttvp72MqcucusRZ200C\/fv3t+bNm9vEiRNrbQszIIAAAggggEBuAgS03NxivdRzzz1ny5cvtwMHDlizZs3s0UcfdWPNtm3bZj\/60Y9s\/fr1Z7VfAU29bF27dq1223SJUzcH8HDaWO9+GocAAgggUAICBLQS2InVbcLx48fPutMyyuYS0KLosSwCCCCAAAKZCxDQMrcKfk4FtH379nEHZ\/CVAAACCCCAQKEFCGiFFi6h9ZeVlVn79u15vVMJ7VM2BQEEEEAgngIEtHjul1i2ioAWy91CoxBAAAEESlCAgFaCO7VQm6QbBPSmAO7gLJQw60UAAQQQQOA\/AgQ0KiFjgaZNm9ratWutV69eGS\/DjAgggAACCCCQvQABLXuzYJfQJc533nkn2O1nwxFAAAEEEKgrAQJaXUmXwOeEVCwlsLvYBAQQQAABjwVCOucG8bL0QtZiSMVSSEfWjQACCCCAQG0CIZ1zCWi1VUMt\/x5SsUSkYnEEEEAAAQQiCYR0ziWgRSqVsO4oiUjF4ggggAACCEQSIKBF4gtr4ZCKJaw9y9YigAACCMRNIKRzLj1oEasvpGKJSMXiCCCAAAIIRBII6ZxLQItUKlzijMjH4ggggAACCGQsQEDLmIoZQyoW9jYCCCCAAALFFAjpnEsPWsRKC6lYIlKxOAIIIIAAApEEQjrnEtAilQqXOCPysTgCCCCAAAIZCxDQMqZixpCKhb2NAAIIIIBAMQVCOufSgxax0kIqlohULI4AAggggEAkgZDOuQS0SKXCJc6IfCyOAAIIIIBAxgIEtIypmDGkYmFvI4AAAgggUEyBkM65wfSgnTx50lauXGm7d++2rl27Wr9+\/axJkyZn1dm7775ry5Yts4MHD1qHDh2sf\/\/+dt5551VbjyEVSzG\/lHw2AggggAACIZ1zgwloo0ePtiNHjtjw4cNt165dtmLFCtuwYUOlkHb48GHr06ePjRgxwr7yla\/Yxo0bbe3atbZp0yZr2LBh2m9GSMXCoQEBBBBAAIFiCoR0zg0ioB04cMCGDh1qW7ZsscaNG7vaGjlypPXt29cGDhyYrLXVq1e70DZ37lz3t9OnT1uPHj1s8eLFrjct3RRSsRTzS8lnI4AAAgggENI5N4iAtmrVKtu6datNnz69Uhjbv3+\/jR07Nvm3jz76yI4dO2af\/\/zn3d8Swe7FF1+0Ro0aVRvQavrKjBo1ysrLy\/lWIYAAAggggEAtArNnz7ZZs2bVONfevXuDcAwioM2bN89OnTpVKSht27bNFi1alOwtq7q3X3jhBbv\/\/vtt4sSJdsMNN1RbDCGl+SC+EWwkAggggEBsBUI65wYR0JYuXWoaX5bak6XxZbqkOWfOnEqFqCD3wAMPmHrXpk6dahdffHGNhRpSscT2G0vDEEAAAQSCEAjpnBtEQFNv2ZIlS0xdp4lp\/vz57qaBcePGVSrqKVOmuP9\/7733VntjQOoCIRVLEN9+NhIBBBBAILYCIZ1zgwhoR48edTcE6FJn586draKiwoYMGWILFixwg\/937txpHTt2tO3bt9sjjzxia9assQYNGrgCPXPmjH388cfJ\/1+1akMqlth+Y2kYAggggEAQAiGdc4MIaKpaPVpj4cKFtmfPHtdzNn78eHcH5\/Hjx92dmk8\/\/bR79Mby5cutadOmrtB1uVOXRmfMmFHtOLSQiiWIbz8biQACCCAQW4GQzrnBBLREtSmQpXtAba7VGFKx5GrEcggggAACCORDIKRzbnABLR8FkrqOkIol33asDwEEEEAAgWwEQjrnEtCyqYw084ZULBGpWBwBBBBAAIFIAiGdcwlokUrFLKRiiUjF4ggggAACCEQSCOmcS0CLVCoEtIh8LI4AAggggEDGAgS0jKmYMaRiYW8jgAACCCBQTIGQzrn0oEWstJCKJSIViyOAAAIIIBBJIKRzLgEtUqlwiTMiH4sjgAACCCCQsQABLWMqZgypWNjbCCCAAAIIFFMgpHMuPWgRKy2kYolIxeIIIIAAAghEEgjpnEtAi1QqXOKMyMfiCCCAAAIIZCxAQMuYihlDKhb2NgIIIIAAAsUUCOmcSw9axEoLqVgiUrE4AggggAACkQRCOucS0CKVCpc4I\/KxOAIIIIAAAhkLENAypmLGkIqFvY0AAggggEAxBUI659KDFrHSQiqWiFQsjgACCCCAQCSBkM65BLRIpcIlzoh8LI4AAggggEDGAgS0jKmYMaRiYW8jgAACCCBQTIGQzrn0oEWstJCKJSIViyOAAAIIIBBJIKRzLgEtUqn4f4lz1qxZVl5eHlGhOIv73HaJ0f7i1A32xXNPfDK1X7x94Ls9Aa14tePdJ\/teLD633+e2q9Bpf\/G+7tgXz57axz6KgO\/f3Wy2nR60Klrvv\/++PfPMM1ZRUWHXX3+9devWzRo0aFCtqe\/F4nP7fW47J6lsDlP5n5fayb9pNmvEPxut\/M6LfX49C7k2AlqK7pEjR+ymm26ynj17Wu\/evW3FihVWv359mz59OgGtkFWY47o50OQIl6fFfPb3ue2E+zwVcITV+Fw\/Pre9FGo\/m7IjoKVorVy50rZu3WozZ850fz1x4oRde+219tRTT1nbtm3TulLs2ZRbfufFPr+e2a7NZ3+f214KJyn8s\/225W9+7PNnWeg1EdBShMePH2+dO3e2IUOGJP86YcIEd6lTQS3dNGzYMNu+fXuh9xPrRwABBBBAIHiBq666ypYsWRKEAwEtZTcPHz7c7rjjDuvevXvyr7rjpUWLFqZ\/Y0IAAQQQQAABBOpCgICWonzPPfe43rPUgDZp0iS78MIL7ZZbbqmL\/cFnIIAAAggggAACRkBLKYI5c+ZYy5YtK4WxW2+91fWeXXfddZQLAggggAACCCBQJwIEtBTmnTt3msahPfnkk9a6dWt7+eWXbcyYMbZlyxZr3LhxnewQPgQBBBBAAAEEECCgVamB5557zpYvX24HDhywZs2a2aOPPmqXX345lYIAAggggAACCNSZAAGtGurjx49bkyZN6mxH8EEIIIAAAggggEBCgIBGLSCAAAIIIIAAAjETIKDluEPeeOMNW7dunZ05c8YGDBjg3qvo07R48WL75z\/\/6W6K+O53v+tT011b\/\/SnP9maNWtMPZ1f+tKXbNCgQV71eL7wwgv2hz\/8wf7973+bnuuj9terV8+r\/fDWW2+5\/XDzzTd70+53333Xli1bZo0aNXJtPn36tH3xi1+0r3\/9695sg15Hp7ecaBhGly5d7MYbb\/Si9vXIopMnTyZfnSd7\/ae8vNyL9qtANE5Z39sPPvjAPTNTb56p6VWAcSuqTZs2ufafOnXKBg8ebF\/+8pfj1sRK7VF7\/\/znP6etE9WSHi6\/e\/du69q1q\/Xr18+bOsoUnYCWqVTKfCoIhRr959JLL7VHHnnE3emZ+oDbHFZbp4u89NJLtmfPHlu\/fr17U4JP07Zt2+zOO+80PQJFz6hbtGiRffrTn7af\/vSnXmzG6tWr7Ze\/\/KVr\/7Fjx2zatGnuQO\/To1wUjL\/xjW\/YBRdcYE888YQX7olgr++rnneoSQf5z33uc3bllVd6sQ2HDx+2\/v37u1pRsJw\/f75ddNFFrpbiPml8rya9Pk+TTry\/\/\/3vTX\/\/xCc+Effm2759++w73\/mO\/eAHP3DHHd31r4A8duzY2LddDVy1apXNnTvXfvzjH9vBgwdNgVnHzNTHSsVtQ958803Tj6q77rrL\/RjUuPDENHr0aNPrGXXu3bVrl\/vRsmHDhpIKaQS0HCpy3Lhx1qlTp+QJ9a9\/\/at973vfsxdffDGHtRVvkb\/85S\/uPaMKOD5N999\/v11yySX27W9\/2zV7\/\/799q1vfcsd8H2YFNDatGmTPDDq\/+skpZOtL9N9991nCgvqQV6wYIEvzXaB4I9\/\/KM99NBD3rQ5taHq+dYPK4VMTeoFVy+DAr5Pk3qO1Wt877332le\/+lUvmv6zn\/3M\/RDUsV7Tjh07XB397ne\/86L9N9xwg02ePDn5Y2Tt2rW2dOlS+\/Wvfx379usKleo+0Vup3uOhQ4dWesLCyJEjrW\/fvjZw4MDYb0+mDSSgZSqVMl+3bt3cqybUe5aY9LeNGzda8+bNc1hjcRZRQJsxY4ZXPSCS0iWeT37yk+7yrCb9MlQ41rb4NulENXXqVHfJTcHfh0kHdJkrKKsn8LHHHvOh2a6N+pWt+tFlQQV7XV726RE6Cgfq\/VPP5WuvveaOQe3atfPGP9FQ1ZDee6xeKF8m\/ZDV0BZ9XzUcQdvw6quvJt\/dHPftUMiR+bnnnuua+vbbb7uQ7MMPW7VdHSGJYSA6\/mhb1MGQmPRDV99pX3o0M6kXAlomSinz6ISqx26oSzX1wK40rxOWT4\/k8DWgJXaH9sXjjz9uCxcudL8E1avmy\/TRRx\/Z3Xff7eqorKzMtT+1+z6u26H26peqDpAffvihO0D6FNB0WUch7ZprrnEBX2MB9exDXx5E3adPHxfONOmkpR5B9ST7dHlcl\/W1Hb\/5zW\/cW1p8mXQ5XJc4\/\/a3v1mrVq1cD7J6zz772c96sQmqkd69e7t60dg\/PUJKx071BMa9Y6FqQJs3b54bR6fxi4lJQ18UonUZt1QmAlqWe1KXdDTu4JVXXqkU0NR9rJOVxoX4Mvkc0NQLoqCgHgSdYBO9ab7YK1y+99579s4777gHI+sAlHqwieN2yFw\/RH7xi1+4HyIaH6Ka18FS3wsfBkvre9u0aVM3REGTDurqefXhMo\/aq4HQPXv2tAceeMC1X2OJdIODtiNx40Mcaye1TRrYvXnzZq96z9R+DUGQs8Y8qYaeffZZ12OjMZgNGzaMO7trq8a76r81dktjSDUGWZfI497+qgFNP2gVkFOPmbqCpV40n3playsaAlptQmn+Xb9E1Ft22WWXuX\/VLyvdSbh9+3avgoKvAe3EiROmV3BpP\/h0912ilHRwTA2UFRUV9s1vftPVT5wPlOp50klKJydN6gn5xz\/+4XoQ9Gtcd1L5NmkMlwKPLl3F2T7hmm6cjYZXaN+oJ9aHSXcPajt86bWU6ccff+yO9xorevHFFztm9eCo5vVg8w4dOvhA79qo3nsNETl06JDrTVNIi\/tUNaApKGuY0ezZs5NN17FJx1ZfhopkYk5Ay0SpyjwaFK07OfVrRJMus+nSj2\/JXbeMz5w507sxaHL\/17\/+lexF0D5Ql73uDov7oyrUzl69etn3v\/99F8rU86THPqimfLvJRKFGvU8+3STw4IMPujsGR40a5QKZDvDaDl\/uRNXJ9Fe\/+pW7jKOxRKoZDbTXeBwfAqZOoBr3V\/WOvBwOw3W+iO6e1c0Y6kXWcUZjt\/S\/U8d11XmjsvjAKVOm2Gc+8xl3mVbHSt2Nqo4FHy6PK6Cpxz5xB\/DRo0fdDQHqvdfjTvQjV09R0LHIp6tYte0+AlptQmn+Xb+cdEDXwVI9CLpcooGjGpfg06RfITrQ+3YXp36B6y6eRE+O9oe6u3Wy8mE8iA40ukSl2tFlKT2LSI9s8Wn8oupcJ1kdIH0KaOox03f3+eefdzWjg\/tPfvITO++887z56mocncYO6Xij74DuzLviiiu8aL\/GO8lbl6J8mzRI\/Yc\/\/GHye6ueNIUCDW\/xYdJwCv0g0WV+Tbr6oAH1cb80rkf69OjRw92xmTpWTp0i+h7o7k4Ffw11KaU7OLWPCGgRvlnq9lY48OkusAiby6J5FtCBR70ecT9A5nmzY7G6xENSfXj+VjowjWHUsYfX0dV9OWlIi+rHV3vVjcaLJnqj6l4wv59Yyq9lJKDlt1ZYGwIIIIAAAgggEFmAgBaZkBUggAACCCCAQKEEdEn89ddfrzTuuFCfFaf1EtDitDdoCwIIIIAAAiUmoMuqeiSPLq1qzLYG+KdONb3bWmN2NU5ad75PmDChxGRq3hwCWlC7m41FAAEEEECgbgU0Zk\/vydTz73QjV+rT\/mt6t7Vu5NG7Q3Vzhu5eJqDV7X7j0xBAAAEEEEAgAAGFLN3goteVJabq3m2tMHfnnXe6Z7Xp5fR67pnevuLTHddRdyk9aFEFWR4BBBBAAAEEahVIF9Cqe7f1unXr7Le\/\/a27y10PxdZjNm677Tb3mrBQJgJaKHua7UQAAQQQQKAAAn\/\/+99t8eLF7uG3iVe+vfXWW+55d6lP9q8a0DJ9t7WeX8glzgLsOFaJAAIIIIAAAqUroKClcWIKZw8\/\/LB736cuY+oVcFdeeWVyw6sGtFJ6t3Uh9i49aIVQZZ0IIIAAAggEJKAbATSgX69h0iMxdNem7thMndJd4iyVd1sXYlcT0AqhyjoRQAABBBAITECPyxg0aJB737DeTV31\/bB6LZzewHP77bcnZUrl3daF2NUEtEKosk4EEEAAAQQCEtDjMvT4jBkzZtizzz5r7733nntHdeqr7ObPn296VVZ5eXlSplTebV2IXU1AK4Qq60QAAQQQQCAQAd0kMGLECJs7d661b9\/ebbWCmgb3T5o0KSMF3m19NhMBLaPSYSYEEEAAAQQQSCegwf6HDh2yc889t9I\/f\/DBB+7BtEy5CRDQcnNjKQQQQAABBBBAoGACBLSC0bJiBBBAAAEEEEAgNwECWm5uLIUAAggggAACCBRMgIBWMFpWjAACCCCAAAII5CZAQMvNjaUQQAABBBBAAIGCCRDQCkbLihFAAAEEEEAAgdwECGi5ubEUAggggAACCCBQMAECWsFoWTECCCCAAAIIIJCbAAEtNzeWQgABBBBAAAEECiZAQCsYLStGAAEEEEAAAQRyEyCg5ebGUggggAACCCCAQMEECGgFo2XFCCCAAAIIIIBAbgIEtNzcWAoBBBBAAAEEECiYAAGtYLSsGAEEEEAAAQQQyE2AgJabG0shgAACCCCAAAIFEyCgFYyWFSOAAAIIIIAAArkJENByc2MpBBBAAAEEEECgYAIEtIi0\/\/M\/\/2Pbt2+31q1bR1wTiyOAAAIIIIBATQJXXXWVLVmyJAikYALayZMnbeXKlbZ7927r2rWr9evXz5o0aXLWTtZ8K1assNcotwWXAAAgAElEQVRee83OP\/98u\/nmm61NmzbVFsM555xj5eXlNnHiRC8L5gtf+ILt3buXthdBwGd7cfncfp\/b7rs97S\/CwSblI6n94vpn8+nBBLTRo0fbkSNHbPjw4bZr1y4XwjZs2HBWSLv99tutRYsWNmjQIFuzZo29+uqr9vzzz1ujRo3SuhLQsim3\/M7LgSa\/ntmuzWd\/n9tOwMm2UvM\/v8\/143PbS6H2s6nGIALagQMHbOjQobZlyxZr3Lix8xk5cqT17dvXBg4cmPR644037O6777Z169a5+c6cOWPLli1z8zRr1iyta1lZmbVv3942bdqUjXts5vX5y+pz20vhQOOzv89tp3aKf\/j0uX58bnsp1H421RtEQFu1apVt3brVpk+fnrRZvXq17d+\/38aOHZv824IFC+z06dM2bNgw27Fjh7Vs2dI6d+5co2eHDh1s3759pv9ON40aNcpdAo3rNHv2bFMbfZx8bru8aX\/xqg774tlT+9jXJKDv5qxZs2pE8nVYTrZ7PoiANm\/ePDt16lSloLRt2zZbtGiRzZ07N2n24IMPujFqmrdHjx6meS644AKbNm2a1a9fP61tIqCpt40JAQQQQAABBAon4HsPYDYyQQS0pUuX2uHDhysFtI0bN5p60ebMmZP0evjhh11Pm8ac1atXz13i\/NrXvuZ63jp16lRtQNM4tMmTJ1uvXr2ysWdeBBBAAAEEEMhCgICWBZYPs6onTLflqus0Mc2fP9\/dNDBu3Ljk3xTk3nzzTZs0aVLybxMmTHB3faaOVUvdZvWg9e\/f35o3b+7tnZw+7EPaiAACCCCAAAGtxGrg6NGj7oYAXerUmLKKigobMmSIacyZAtbOnTutY8eO9uGHH7qwpRsDLrvsMjffjTfe6P6\/bgRINxHQSqxY2BwEEEAAgdgKENBiu2tyb5gerbFw4ULbs2eP6zkbP3686xU7fvy4G2\/29NNPW7t27UyXPh966CH3QRqLph60AQMGVPvBCmhTpkxxvXO+3smZuypLIoAAAgggUHcCBLS6s67zT1IgS\/eA2tSGaOzZsWPHqn20Ruq8CmibN292D7XlRoE63518IAIIIIBAQAIEtIB2dtRNJaBFFWR5BBBAAAEEMhMgoGXmxFxmbgybetB+\/vOfu8dxJO7kTLz6iTs7KRMEEEAAAQTyI0BAy49jEGtJBLTExr7yyivufyqsaUq8uUCBjbAWREmwkQgggAACBRIgoBUIthRXWzWgVd3G1MCm\/53aw0ZgK8WKYJsQQAABBAolQEArlGwJrre2gJZJYFNQu\/baa+lhK8H6YJMQQAABBPInQEDLn2XJrynbgJYusKln7eWXX7ZEDxuBreTLhg1EAAEEEMhBgICWA1qoi0QNaJkEttR50l0WVe9b1YnLp6FWJNuNAAIIlK4AAa10923etyzfAa26S6KJvyfGtKXOp9632pbTv9cW2mr79yh46UJklPWlW7aQ7c93W1kfAggggED2AgS07M2CXaKsrMxuuukmGzNmTOwN0oW71EbX9u9RNjBdiIyyvnTLFrL9mQTcfGxPIUMmITkfe4h1IIBAMQUIaMXU9+yzL730CuvT5wbbteslu\/rqq6272q8erauvrnlLurs53fT\/\/1eVRVLm8YylJJtb6AAotEJ+BiE53mVZyHAe7y2ndYUUqIsfZoVsf9V1\/9\/\/+3\/tnXfeqcuPLNpn1TvD+4ki4ZeVDbM2bR63HTt0mXGLmW1OWd9\/nolW+3Si9lmYwxuBtIGbsJ2X\/VfIAJuXBrISBBAoqMCnPvUp9yrGECYCWsS9\/H\/KLrF\/\/G9n22EbzlpT9+5dslr71VdnN39WK2fmkhLIV2+Yen2LMeWr\/Vm1\/b8Pkc5qmZjPnOlPwJhvBs1DIGMBAlrGVMxYVnaJ3XTTMBtjR81iNg4tXW9DlB6IopxU81xiUbY\/z02JtDrfL4f53v5IOy+PC5fa5as80rCqagR8\/+4xBo3SzligkHdxKkzk+67NKF\/OUjgZRNn+jIuCGRFAAAEECiJAQCsIa2muNF8BLTWMVX1obVU5nntWmrXEViGAAAII1CxAQKNCMhbIJaDVFsZ47VPG\/MyIAAIIIBCQAAEtoJ0ddVNrC2ipYWzatGnu43SZLXGpjTAWdQ+wPAIIIIBAKAIEtFD2dB62M11A+\/nPf27pwtjEiRPz8ImsAgEEEEAAgTAFCGhh7vectrpqQBs8eLD169fPCGM5cbIQAggggAAC1QoQ0CiOjAVSA5ouZ86ePds2bdqU8fLMiAACCCCAAAKZCRDQMnPyaq7333\/fnnnmGauoqLDrr7\/eunXrZg0aNKhxGxYvXmzXXHONtWvXrtr5UgOaes8mT55c60vJvYKjsQgggAACCMREgIAWkx2Rr2YcOXLEvdC8Z8+e1rt3b1uxYoXVr1\/fpk+fXu1HrF692u677z578sknrUuX6p\/wnwho9J7la2+xHgQQQAABBNILENBKrDJWrlxpW7dutZkzZ7otO3HihOnuyaeeesratm171tbu3bvXhg8fbhdddJELaZ07d66xB+306dN28OBBO+ecc6xJkyaV5h01apSVl5eXmCibgwACCCCAQP4FNExo1qxZNa5Y5+gQpiDexTl+\/HgXsoYMGZLcpxMmTHCXOqs+9PXo0aM2cOBAGzdunK1bt85uu+0269SpU40BbfPmzXb++ecb750P4SvDNiKAAAIIFEuAHrRiyRfoc9Ubdscdd1j37t2Tn6CE3qJFC9dTlpgUsNTjdfnll9tdd93l\/veIESNqDWhTpkzh5oAC7TtWiwACCCCAQEKAgFZitXDPPfe43rPUgDZp0iS78MIL7ZZbbklu7Zw5c2z\/\/v3uGWb16tVLBrSOHTtaw4YN06poDBoBrcQKhs1BAAEEEIilAAEtlrsl90YpeLVs2bJSGLv11ltd79l1112XXLGeX6abBxo1auT+puvcrVq1cj1ojz32GAEt913AkggggAACCEQWIKBFJozXCnbu3Gkah6Y7Mlu3bm16GfmYMWNsy5Yt7oaBt99+O+1lTF3m1CXO2m4S0M0BPF4jXvuc1iCAAAIIlJ4AAa309qk999xztnz5cjtw4IA1a9bMHn30UTfWbNu2bfajH\/3I1q9ff9ZWK6Cpl61r167ViugSJwGtBAuGTUIAAQQQiJ0AAS12uyR\/DTp+\/PhZj8KIsnYFtH379nEHZxRElkUAAQQQQCADAQJaBkjM8h8BAhqVgAACCCCAQN0IENDqxrkkPqWsrMzat2\/P+zdLYm+yEQgggAACcRYgoMV578SsbQS0mO0QmoMAAgggULICBLSS3bX53zDdIKBXOU2cODH\/K2eNCCCAAAIIIJAUIKBRDBkLqAdt8eLF1qtXr4yXYUYEEEAAAQQQyF6AgJa9WbBLhFQswe5kNhwBBBBAIBYCIZ1zg3hZeiGrKqRiKaQj60YAAQQQQKA2gZDOuQS02qqhln8PqVgiUrE4AggggAACkQRCOucS0CKVillIxRKRisURQAABBBCIJBDSOZeAFqlUCGgR+VgcAQQQQACBjAUIaBlTMWNIxcLeRgABBBBAoJgCIZ1z6UGLWGkhFUtEKhZHAAEEEEAgkkBI51wCWqRS4RJnRD4WRwABBBBAIGMBAlrGVMwYUrGwtxFAAAEEECimQEjnXHrQIlZaSMUSkYrFEUAAAQQQiCQQ0jmXgBapVLjEGZGPxRFAAAEEEMhYgICWMRUzhlQs7G0EEEAAAQSKKRDSOZcetIiVFlKxRKRicQQQQAABBCIJhHTOJaBFKhUucUbkY3EEEEAAAQQyFiCgZUzFjCEVC3sbAQQQQACBYgqEdM4Npgft5MmTtnLlStu9e7d17drV+vXrZ02aNDmrzt59911btmyZHTx40Dp06GD9+\/e38847r9p6DKlYivml5LMRQAABBBAI6ZwbTEAbPXq0HTlyxIYPH267du2yFStW2IYNGyqFtMOHD1ufPn1sxIgR9pWvfMU2btxoa9eutU2bNlnDhg3TfjNCKhYODQgggAACCBRTIKRzbhAB7cCBAzZ06FDbsmWLNW7c2NXWyJEjrW\/fvjZw4MBkra1evdqFtrlz57q\/nT592nr06GGLFy92vWnpJhVLTdOoUaOsvLy8mPXMZyOAAAIIIOCFwOzZs23WrFk1tnXv3r1ebEvURgYR0FatWmVbt2616dOnVwpj+\/fvt7Fjxyb\/9tFHH9mxY8fs85\/\/vPtbIti9+OKL1qhRI3rQolYbyyOAAAIIIBBBgB60CHhxXHTevHl26tSpSj1Z27Zts0WLFiV7y6q2+4UXXrD777\/fJk6caDfccEO1mxVSscRx39ImBBBAAIFwBEI65wbRg7Z06VLT+LLUS40aX6ZLmnPmzKlU2QpyDzzwgKl3berUqXbxxRfXWPkhFUs4hwC2FAEEEEAgjgIhnXODCGjqLVuyZInp2nZimj9\/vrtpYNy4cZVqcMqUKe7\/33vvvdXeGJC6QEjFEscvK21CAAEEEAhHIKRzbhAB7ejRo+6GAF3q7Ny5s1VUVNiQIUNswYIFbvD\/zp07rWPHjrZ9+3Z75JFHbM2aNdagQQNX8WfOnLGPP\/44+f+rfg1CKpZwDgFsKQIIIIBAHAVCOucGEdBUZHq0xsKFC23Pnj2u52z8+PHuDs7jx4+7OzWffvpp9+iN5cuXW9OmTV1d6nKnLo3OmDGj2nFoIRVLHL+stAkBBBBAIByBkM65wQS0RPkqkKV7QG2u5R1SseRqxHIIIIAAAgjkQyCkc25wAS0fBZK6jpCKJd92rA8BBBBAAIFsBEI65xLQsqmMNPOGVCwRqVgcAQQQQACBSAIhnXMJaJFKxSykYolIxeIIIIAAAghEEgjpnEtAi1QqBLSIfCyOAAIIIIBAxgIEtIypmDGkYmFvI4AAAgggUEyBkM659KBFrLSQiiUiFYsjgAACCCAQSSCkcy4BLVKpcIkzIh+LI4AAAgggkLEAAS1jKmYMqVjY2wgggAACCBRTIKRzLj1oESstpGKJSMXiCCCAAAIIRBII6ZxLQItUKlzijMjH4ggggAACCGQsQEDLmIoZQyoW9jYCCCCAAALFFAjpnEsPWsRKC6lYIlKxOAIIIIAAApEEQjrnEtAilQqXOCPysTgCCCCAAAIZCxDQMqZixpCKhb2NAAIIIIBAMQVCOufSgxax0kIqlohULI4AAggggEAkgZDOuQS0SKXCJc6IfCyOAAIIIIBAxgIEtIypmDGkYmFvI4AAAgggUEyBkM659KBFrDTfi2XWrFlWXl4eUaE4i\/vcdonR\/uLUDfbFc098MrVfvH3gu73v59xs9jwBLRutNPP6Xiw+t9\/ntquUaH\/EL1+ExbGPgJeHRfHPA2KOq8A+R7giLEZAq4L+\/vvv2zPPPGMVFRV2\/fXXW7du3axBgwbV7hqKvQhV+9+PxL549r4HTGqH2oki4HP9+Nx234872dYcAS1F7MiRI3bTTTdZz549rXfv3rZixQqrX7++TZ8+nYCWbWXVwfwcaOoAuYaP8Nnf57aXwkkK\/+J9d7Evnn22n0xASxFbuXKlbd261WbOnOn+euLECbv22mvtqaeesrZt26a1pdizLbn8zY99\/ixzWZPP\/j63nYCWS7Xmdxmf68fntpdC7WdTiQS0FK3x48db586dbciQIcm\/TpgwwV3qVFBLNw0bNsy2b9+ejTnzIoAAAggggEAOAldddZUtWbIkhyX9W4SAlrLPhg8fbnfccYd17949+Vfd8dKiRQvTvzEhgAACCCCAAAJ1IUBAS1G+5557XO9ZakCbNGmSXXjhhXbLLbfUxf7gMxBAAAEEEEAAASOgpRTBnDlzrGXLlpXC2K233up6z6677jrKBQEEEEAAAQQQqBMBAloK886dO03j0J588klr3bq1vfzyyzZmzBjbsmWLNW7cuE52CB+CAAIIIIAAAggQ0KrUwHPPPWfLly+3AwcOWLNmzezRRx+1yy+\/nEpBAAEEEEAAAQTqTICAVg318ePHrUmTJnW2I\/ggBBBAAAEEEEAgIUBAoxYQQAABBBBAAIGYCRDQctwhb7zxhq1bt87OnDljAwYMcO9V9GlavHix\/fOf\/3Q3RXz3u9\/1qemurX\/6059szZo1pp7OL33pSzZo0CCvejxfeOEF+8Mf\/mD\/\/ve\/Tc\/1Ufvr1avn1X5466233H64+eabvWn3u+++a8uWLbNGjRq5Np8+fdq++MUv2te\/\/nVvtkGvo9NbTjQMo0uXLnbjjTd6Uft6ZNHJkyeTr86Tvf5TXl7uRftVIBqnrO\/tBx984J6ZqTfP1PQqwLgV1aZNm1z7T506ZYMHD7Yvf\/nLcWtipfaovX\/+85\/T1olqSQ+X3717t3Xt2tX69evnTR1lik5Ay1QqZT4VhEKN\/nPppZfaI4884u70TH3AbQ6rrdNFXnrpJduzZ4+tX7\/evSnBp2nbtm125513mh6BomfULVq0yD796U\/bT3\/6Uy82Y\/Xq1fbLX\/7Stf\/YsWM2bdo0d6D36VEuCsbf+MY37IILLrAnnnjCC\/dEsNf3Vc871KSD\/Oc+9zm78sorvdiGw4cPW\/\/+\/V2tKFjOnz\/fLrroIldLcZ80vleTXp+nSSfe3\/\/+96a\/f+ITn4h7823fvn32ne98x37wgx+4447u+ldAHjt2bOzbrgauWrXK5s6daz\/+8Y\/t4MGDpsCsY2bqY6XitiFvvvmm6UfVXXfd5X4Malx4Yho9erTp9Yw69+7atcv9aNmwYUNJhTQCWg4VOW7cOOvUqVPyhPrXv\/7Vvve979mLL76Yw9qKt8hf\/vIX955RBRyfpvvvv98uueQS+\/a3v+2avX\/\/fvvWt77lDvg+TApobdq0SR4Y9f91ktLJ1pfpvvvuM4UF9SAvWLDAl2a7QPDHP\/7RHnroIW\/anNpQ9Xzrh5VCpib1gquXQQHfp0k9x+o1vvfee+2rX\/2qF03\/2c9+5n4I6livaceOHa6Ofve733nR\/htuuMEmT56c\/DGydu1aW7p0qf3617+Offt1hUp1n+itVO\/x0KFDKz1hYeTIkda3b18bOHBg7Lcn0wYS0DKVSpmvW7du7lUT6j1LTPrbxo0brXnz5jmssTiLKKDNmDHDqx4QSekSzyc\/+Ul3eVaTfhkqHGtbfJt0opo6daq75Kbg78OkA7rMFZTVE\/jYY4\/50GzXRv3KVv3osqCCvS4v+\/QIHYUD9f6p5\/K1115zx6B27dp5459oqGpI7z1WL5Qvk37IamiLvq8ajqBtePXVV5Pvbo77dijkyPzcc891TX377bddSPbhh63aro6QxDAQHX+0LepgSEz6oavvtC89mpnUCwEtE6WUeXRC1WM31KWaemBXmtcJy6dHcvga0BK7Q\/vi8ccft4ULF7pfgupV82X66KOP7O6773Z1VFZW5tqf2n0f1+1Qe\/VLVQfIDz\/80B0gfQpouqyjkHbNNde4gK+xgHr2oS8Pou7Tp48LZ5p00lKPoHqSfbo8rsv62o7f\/OY37i0tvky6HK5LnH\/729+sVatWrgdZvWef\/exnvdgE1Ujv3r1dvWjsnx4hpWOnegLj3rFQNaDNmzfPjaPT+MXEpKEvCtG6jFsqEwEtyz2pSzoad\/DKK69UCmjqPtbJSuNCfJl8DmjqBVFQUA+CTrCJ3jRf7BUu33vvPXvnnXfcg5F1AEo92MRxO2SuHyK\/+MUv3A8RjQ9Rzetgqe+FD4Ol9b1t2rSpG6KgSQd19bz6cJlH7dVA6J49e9oDDzzg2q+xRLrBQduRuPEhjrWT2iYN7N68ebNXvWdqv4YgyFljnlRDzz77rOux0RjMhg0bxp3dtVXjXfXfGrulMaQag6xL5HFvf9WAph+0Csipx0xdwVIvmk+9srUVDQGtNqE0\/65fIuotu+yyy9y\/6peV7iTcvn27V0HB14B24sQJ0yu4tB98uvsuUUo6OKYGyoqKCvvmN7\/p6ifOB0r1POkkpZOTJvWE\/OMf\/3A9CPo1rjupfJs0hkuBR5eu4myfcE03zkbDK7Rv1BPrw6S7B7UdvvRayvTjjz92x3uNFb344osds3pwVPN6sHmHDh18oHdtVO+9hogcOnTI9aYppMV9qhrQFJQ1zGj27NnJpuvYpGOrL0NFMjEnoGWiVGUeDYrWnZz6NaJJl9l06ce35K5bxmfOnOndGDS5\/+tf\/0r2ImgfqMted4fF\/VEVamevXr3s+9\/\/vgtl6nnSYx9UU77dZKJQo94nn24SePDBB90dg6NGjXKBTAd4bYcvd6LqZPqrX\/3KXcbRWCLVjAbaazyODwFTJ1CN+6t6R14Oh+E6X0R3z+pmDPUi6zijsVv636njuuq8UVl84JQpU+wzn\/mMu0yrY6XuRlXHgg+XxxXQ1GOfuAP46NGj7oYA9d7rcSf6kaunKOhY5NNVrNp2HwGtNqE0\/65fTjqg62CpHgRdLtHAUY1L8GnSrxAd6H27i1O\/wHUXT6InR\/tD3d06WfkwHkQHGl2iUu3ospSeRaRHtvg0flF1rpOsDpA+BTT1mOm7+\/zzz7ua0cH9Jz\/5iZ133nnefHU1jk5jh3S80XdAd+ZdccUVXrRf453krUtRvk0apP7DH\/4w+b1VT5pCgYa3+DBpOIV+kOgyvyZdfdCA+rhfGtcjfXr06OHu2EwdK6dOEX0PdHengr+GupTSHZzaRwS0CN8sdXsrHPh0F1iEzWXRPAvowKNej7gfIPO82bFYXeIhqT48fysdmMYw6tjD6+jqvpw0pEX146u96kbjRRO9UXUvmN9PLOXXMhLQ8lsrrA0BBBBAAAEEEIgsQECLTMgKEEAAAQQQQACB\/AoQ0PLrydoQQAABBBBAII8CGrP4+uuvV7oxLI+rj+2qCGix3TU0DAEEEEAAAf8FNO5Nz0zU2DfdVKc7MFMn3Um9bt06d1f7gAED3HMhE5NuqtKNbHo00YQJE\/zHyGILCGhZYDErAggggAACCGQnoJsq9CJzPaBYd9qnvo5Jj6zSXez6jx48rvfM6mHAukNWd1rr5e66e1aPlyGgZefO3AgggAACCCCAQK0CClm6A1nvk01MerCsetUSz2PT40z0zlmFuTvvvNM9TLdFixbuwbR6PZ5Pj8SpFaSWGehBiyrI8ggggAACCCBQq0C6gKY3YSh8qfcsMelvuuT529\/+1j2GSG8t0XPQbrvtNvce11AmAlooe5rtRAABBBBAoAACf\/\/7323x4sXu7QSJd\/K+9dZb7oHEqa9eqhrQ1JumB3TrobOpzxPVGxr0OsXEw7v1gGkucRZgx7FKBBBAAAEEEChdAQUtjRNTOHv44YfdC9l1GVPv6L3yyiuTG141oOmmgC5duri3G6QGNL2dQTcVlNJrm3LZ+\/Sg5aLGMggggAACCCCQFNCNABrQr\/dk6pEYClgaW5Y6pbvEqbFn6i3Ty+g16U0Nekfo9u3b3Z2bIU8EtJD3PtuOAAIIIIBAngT0uIxBgwZZr169bM6cOe5VdqmT3turVyTefvvtyT\/rb7qTc9q0ae5vjz\/+uLvkqeVDnwhooVcA248AAggggEBEAYUsPT5jxowZ9uyzz9p7771nU6dOrfSu4fnz57sesvLy8uSn6RlpTzzxhK1fv969iF69blquVatWEVvk\/+IENP\/3IVuAAAIIIIBA0QR0k8CIESNs7ty51r59e9cOBTUN7p80aVJG7VLPmsJa6li0jBYs4ZkIaCW8c9k0BBBAAAEECi2gwf6HDh2yc889t9JHffDBB+7BtEy5CRDQcnNjKQQQQAABBBBAoGACBLSC0bJiBBBAAAEEEEAgNwECWm5uLIUAAggggAACCBRMgIBWMFpWjAACCCCAAAII5CZAQMvNjaUQQAABBBBAAIGCCRDQCkbLihFAAAEEEEAAgdwECGi5ubEUAggggAACCCBQMAECWsFoWTECCCCAAAIIIJCbAAEtNzeWQgABBBBAAAEECiZAQCsYLStGAAEEEEAAAQRyEyCg5ebGUggggAACCCCAQMEECGgFo2XFCCCAAAIIIIBAbgIEtNzcWAoBBBBAAAEEECiYAAGtYLSsGAEEEEAAAQQQyE2AgJabG0shgAACCCCAAAIFEyCgFYyWFSOAAAIIIIAAArkJENByc0suNWzYMNu+fXvEtbA4AggggAACCNQmcNVVV9mSJUtqm60k\/p2AFnE3nnPOOVZeXm4TJ06MuKbiLP6FL3zB9u7dW5wPj\/ipPrddm077IxZAhMWxj4CXh0XxzwNijqvAPke4IiwWTEA7efKkrVy50nbv3m1du3a1fv36WZMmTc4i13wrVqyw1157zc4\/\/3y7+eabrU2bNtXuGgJaEar2vx\/JgaZ49r4HTGqH2oki4HP9+Nx234872dZcMAFt9OjRduTIERs+fLjt2rXLhbANGzacFdJuv\/12a9GihQ0aNMjWrFljr776qj3\/\/PPWqFGjtLYEtGxLLn\/zc6DJn2Uua\/LZ3+e2l8JJCv9cvnH5WQb7\/DjWxVqCCGgHDhywoUOH2pYtW6xx48bOdeTIkda3b18bOHBg0vmNN96wu+++29atW+fmO3PmjC1btszN06xZs7T7o6yszD788ENr27Zt2n8fNWqUuwQa18nnL6vPbeckW9xvBLWDfxQBn+sn7m2fPXu2zZo1q8bd4+uwnGxrLoiAtmrVKtu6datNnz496bN69Wrbv3+\/jR07Nvm3BQsW2OnTp00D\/3fs2GEtW7a0zp0712jaoUMH27dvnwtzPk76MihE+jj53HZ50\/7iVR32xbOn9rGPIhD3gBll26ouG0RAmzdvnp06dapST9a2bdts0aJFNnfu3KTJgw8+6Maoad4ePXqY5rngggts2rRpVr9+\/bTuvge0fBYT60IAAQQQQKCQAgS0QuoWYd1Lly61w4cPVwpoGzduNPWizZkzJ9mihx9+2PW0acxZvXr1XK\/Y1772Ndfz1qlTp2oDmsahTZ482Xr16lWErfAaAfoAABuNSURBVOMjEUAAAQQQCEOAgFZi+1k9YXpuii5rJKb58+e7mwbGjRuX\/JuC3JtvvmmTJk1K\/m3ChAnurs\/UsWqpPOpB69+\/vzVv3tzbR22U2O5mcxBAAAEESlSAgFZiO\/bo0aPuhgBd6tSYsoqKChsyZIhpzJkC1s6dO61jx45usL\/Clm4MuOyyy9x8N954o\/v\/7du3r7YHjYBWYgXD5iCAAAIIxFKAgBbL3RKtUXq0xsKFC23Pnj2u52z8+PGuV+z48eNuvNnTTz9t7dq1M136fOihh9yHaSyaetAGDBhQ7Ycr4E2ZMsX1zm3atClaI1kaAQQQQAABBKoVIKCVcHEokKV7QG3qJmvs2bFjx6p9tEbqvApomzdvdg+19fVOzhLe3WwaAggggEAJCRDQSmhnFnpTCGiFFmb9CCCAAAII\/EeAgEYlZCyQCGg\/\/\/nP3eM4NOluzsS7ObmzM2NKZkQAAQQQQKBGAQIaBZKxQCKgpS7wyiuvJMOa\/jeBLWNOZkQAAQQQQKBaAQIaxZGxQLqAVnVhAlvGnMyIAAIIIIAAAc3MgniTQCFrPZOAVlNgS7wblEuihdxLrBsBBBBAoBQE6EErhb1YR9uQS0CrejlU\/z8xfi01sDF+rY52Ih+DAAIIIOCFAAHNi90Uj0ZGDWjpetcSgU2XRjWlC2rp\/nbttdeehULIi0ed0AoEEEAAgegCBLTohsGsId8BLR1cIqil63lL\/dvLL7981uLplo0a2qIuX1NxpAuZ+S6mQrY\/321lfQgggAAC\/1+AgEY1ZCxQVlZmN910k40ZMybjZYo9Y7rQlk2boi5f02elC5nZtC2TeQvZ\/up6PDNpVzbzFDJkEpKz2RPMiwACdSlAQKtLbc8\/69JLr7A+fW6wXbteclty9dVXW\/c0PVl29dU1b2n37hlLZD5nxqss7oxZbHtxG5rZpxc6AKoVhfwMQnJm+5m5ogsU8odG9NaxhjgKPP744\/bOO+\/EsWl5bxN3cUYkLSsbZv\/7v0vMbLNOm\/\/97y0R11rb4idqm4F\/RwABBBBAoOQEPvWpT7lXMYYwEdAi7uX\/U3aJtWlzuc3fscm6d++S8dquvjrzeWtbafcS64GqbXt9+\/dC9nbVlUVd9KrVxbaUwr6oCyc+A4G4ChDQ4rpnYtiu\/4xBu8XGKJzFJChFuakgbsSlcEIthcs4pbANqu26GF8Xt+8Q7clMoFRqPLOt9XcuxqD5u+\/qvOWFuotTwaRqOIlyl2amj+Woc8BaPpCDZtz2CO1BAAEEiidAQCuevXefnI+AlhrGFMIS7++sGk54zpl35UGDEUAAAQTyKEBAyyNmqa8q24BWWxhTCKPXqNSrhu1DAAEEEMhFgICWi1qgy9QU0AhjgRYFm40AAgggUBABAlpBWEtzpekC2uDBg8+6TEnPWGnuf7YKAQQQQKDuBAhodWft\/SdVDWgKZ\/369bOJEyd6v21sAAIIIIAAAnESIKDFaW\/EvC2pAU2XNGfPnm2bNm2KeatpHgIIIIAAAv4JEND822e1tvj999+3Z555xioqKuz666+3bt26WYMGDWpcbvHixXbNNddYu3btqp0vNaCp92zy5MkM8q91bzADAggggAAC2QsQ0LI3i\/USR44ccS8079mzp\/Xu3dtWrFhh9evXt+nTp1fb7tWrV9t9991nTz75pHXpUv1T\/xMBjd6zWJcAjUMAAQQQKAEBAloJ7MTUTVi5cqVt3brVZs6c6f584sQJ90Txp556ytq2bXvW1u7du9eGDx9uF110kQtpnTt3rrEH7fTp03bw4EE755xzrEmTJpXmHTVqlJWXl5eYKJuDAAIIIIBA\/gU0TGjWrFk1rljn6BCmIN7FOX78eBeyhgwZktynEyZMcJc6qz789ejRozZw4EAbN26crVu3zm677Tbr1KlTjQFt8+bNdv7559uZM2dCqBm2EQEEEEAAgaII0INWFPbCfah6w+644w5Lfam4EnqLFi1cT1liUsBSj9fll19ud911l\/vfI0aMqDWgTZkyhZsDCrf7WDMCCCCAAAJOgIBWYoVwzz33uN6z1IA2adIku\/DCC+2WW25Jbu2cOXNs\/\/79Nm3aNKtXr14yoHXs2NEaNmyYVkVj0AhoJVYwbA4CCCCAQCwFCGix3C25N0rBq2XLlpXC2K233up6z6677rrkivX8Mt080KhRI\/c3Xedu1aqV60F77LHHqg1o\/fv3t+bNm\/Pss9x3EUsigAACCCBQqwABrVYiv2bYuXOnaRya7shs3bq16YXkY8aMsS1btrgbBt5+++20lzF1mVOXOGu7SUA3B\/B4Db9qgtYigAACCPgnQEDzb5\/V2uLnnnvOli9fbgcOHLBmzZrZo48+6saabdu2zX70ox\/Z+vXrz1qHApp62bp27Vrt+nWJk4BWKz8zIIAAAgggEFmAgBaZML4rOH78+FmPwojSWgW0ffv2cQdnFESWRQABBBBAIAMBAloGSMzyHwECGpWAAAIIIIBA3QgQ0OrGuSQ+payszNq3b8\/7N0tib7IRCCCAAAJxFiCgxXnvxKxtBLSY7RCagwACCCBQsgIEtJLdtfnfsKZNm9ratWt5QXr+aVkjAggggAAClQQIaBRExgLqQVu8eDEBLWMxZkQAAQQQQCA3AQJabm5BLhVSsQS5g9loBBBAAIHYCIR0zg3iZemFrKyQiqWQjqwbAQQQQACB2gRCOucS0Gqrhlr+PaRiiUjF4ggggAACCEQSCOmcS0CLVCpmIRVLRCoWRwABBBBAIJJASOdcAlqkUiGgReRjcQQQQAABBDIWIKBlTMWMIRULexsBBBBAAIFiCoR0zqUHLWKlhVQsEalYHAEEEEAAgUgCIZ1zCWiRSoVLnBH5WBwBBBBAAIGMBQhoGVMxY0jFwt5GAAEEEECgmAIhnXPpQYtYaSEVS0QqFkcAAQQQQCCSQEjnXAJapFLhEmdEPhZHAAEEEEAgYwECWsZUzBhSsbC3EUAAAQQQKKZASOdcetAiVlpIxRKRisURQAABBBCIJBDSOZeAFqlUuMQZkY\/FEUAAAQQQyFiAgJYxFTOGVCzsbQQQQAABBIopENI5N5getJMnT9rKlStt9+7d1rVrV+vXr581adLkrDp79913bdmyZXbw4EHr0KGD9e\/f384777xq6zGkYinml5LPRgABBBBAIKRzbjABbfTo0XbkyBEbPny47dq1y1asWGEbNmyoFNIOHz5sffr0sREjRthXvvIV27hxo61du9Y2bdpkDRs2TPvNCKlYODQggAACCCBQTIGQzrlBBLQDBw7Y0KFDbcuWLda4cWNXWyNHjrS+ffvawIEDk7W2evVqF9rmzp3r\/nb69Gnr0aOHLV682PWmpZtULDVNo0aNsvLy8mLWM5+NAAIIIICAFwKzZ8+2WbNm1djWvXv3erEtURsZREBbtWqVbd261aZPn14pjO3fv9\/Gjh2b\/NtHH31kx44ds89\/\/vPub4lg9+KLL1qjRo3oQYtabSyPAAIIIIBABAF60CLgxXHRefPm2alTpyr1ZG3bts0WLVqU7C2r2u4XXnjB7r\/\/fps4caLdcMMN1W5WSMUSx31LmxBAAAEEwhEI6ZwbRA\/a0qVLTePLUi81anyZLmnOmTOnUmUryD3wwAOm3rWpU6faxRdfXGPlh1Qs4RwC2FIEEEAAgTgKhHTODSKgqbdsyZIlpmvbiWn+\/PnupoFx48ZVqsEpU6a4\/3\/vvfdWe2NA6gIhFUscv6y0CQEEEEAgHIGQzrlBBLSjR4+6GwJ0qbNz585WUVFhQ4YMsQULFrjB\/zt37rSOHTva9u3b7ZFHHrE1a9ZYgwYNXMWfOXPGPv744+T\/r\/o1CKlYwjkEsKUIIIAAAnEUCOmcG0RAU5Hp0RoLFy60PXv2uJ6z8ePHuzs4jx8\/7u7UfPrpp92jN5YvX25NmzZ1danLnbo0OmPGjGrHoYVULHH8stImBBBAAIFwBEI65wYT0BLlq0CW7gG1uZZ3SMWSqxHLIYAAAgggkA+BkM65wQW0fBRI6jpCKpZ827E+BBBAAAEEshEI6ZxLQMumMtLMG1KxRKRicQQQQAABBCIJhHTOJaBFKhWzkIolIhWLI4AAAgggEEkgpHMuAS1SqRDQIvKxOAIIIIAAAhkLENAypmLGkIqFvY0AAggggEAxBUI659KDFrHSQiqWiFQsjgACCCCAQCSBkM65BLRIpcIlzoh8LI4AAggggEDGAgS0jKmYMaRiYW8jgAACCCBQTIGQzrn0oEWstJCKJSIViyOAAAIIIBBJIKRzLgEtUqlwiTMiH4sjgAACCCCQsQABLWMqZgypWNjbCCCAAAIIFFMgpHMuPWgRKy2kYolIxeIIIIAAAghEEgjpnEtAi1QqXOKMyMfiCCCAAAIIZCxAQMuYihlDKhb2NgIIIIAAAsUUCOmcSw9axEoLqVgiUrE4AggggAACkQRCOucS0CKVCpc4I\/KxOAIIIIAAAhkLENAypmLGkIqFvY0AAggggEAxBUI659KDFrHSfC+WWbNmWXl5eUSF4izuc9slRvuLUzfYF8898cnUfvH2ge\/2vp9zs9nzBLRstNLM63ux+Nx+n9uuUqL9Eb98ERbHPgJeHhbFPw+IOa4C+xzhirAYAa0K+vvvv2\/PPPOMVVRU2PXXX2\/dunWzBg0aVLtrKPYiVO1\/PxL74tn7HjCpHWonioDP9eNz230\/7mRbcwS0FLEjR47YTTfdZD179rTevXvbihUrrH79+jZ9+nQCWraVVQfzc6CpA+QaPsJnf5\/bXgonKfyL993Fvnj22X4yAS1FbOXKlbZ161abOXOm++uJEyfs2muvtaeeesratm2b1pZiz7bk8jc\/9vmzzGVNPvv73HYCWi7Vmt9lfK4fn9teCrWfTSUS0FK0xo8fb507d7YhQ4Yk\/zphwgR3qVNBLd00bNgw2759ezbmzIsAAggggAACOQhcddVVtmTJkhyW9G8RAlrKPhs+fLjdcccd1r179+RfdcdLixYtTP\/GhAACCCCAAAII1IUAAS1F+Z577nG9Z6kBbdKkSXbhhRfaLbfcUhf7g89AAAEEEEAAAQSMgJZSBHPmzLGWLVtWCmO33nqr6z277rrrKBcEEEAAAQQQQKBOBAhoKcw7d+40jUN78sknrXXr1vbyyy\/bmDFjbMuWLda4ceM62SF8CAIIIIAAAgggQECrUgPPPfecLV++3A4cOGDNmjWzRx991C6\/\/HIqBQEEEEAAAQQQqDMBAlo11MePH7cmTZrU2Y7ggxBAAAEEEEAAgYQAAY1aQAABBBBAAAEEYiZAQMtxh7zxxhu2bt06O3PmjA0YMMC9V9GnafHixfbPf\/7T3RTx3e9+16emu7b+6U9\/sjVr1ph6Or\/0pS\/ZoEGDvOrxfOGFF+wPf\/iD\/fvf\/zY910ftr1evnlf74a233nL74eabb\/am3e+++64tW7bMGjVq5Np8+vRp++IXv2hf\/\/rXvdkGvY5ObznRMIwuXbrYjTfe6EXt65FFJ0+eTL46T\/b6T3l5uRftV4FonLK+tx988IF7ZqbePFPTqwDjVlSbNm1y7T916pQNHjzYvvzlL8etiZXao\/b++c9\/TlsnqiU9XH737t3WtWtX69evnzd1lCk6AS1TqZT5VBAKNfrPpZdeao888oi70zP1Abc5rLZOF3nppZdsz549tn79evemBJ+mbdu22Z133ml6BIqeUbdo0SL79Kc\/bT\/96U+92IzVq1fbL3\/5S9f+Y8eO2bRp09yB3qdHuSgYf+Mb37ALLrjAnnjiCS\/cE8Fe31c971CTDvKf+9zn7Morr\/RiGw4fPmz9+\/d3taJgOX\/+fLvoootcLcV90vheTXp9niadeH\/\/+9+b\/v6JT3wi7s23ffv22Xe+8x37wQ9+4I47uutfAXns2LGxb7sauGrVKps7d679+Mc\/toMHD5oCs46ZqY+VituGvPnmm6YfVXfddZf7Mahx4Ylp9OjRptcz6ty7a9cu96Nlw4YNJRXSCGg5VOS4ceOsU6dOyRPqX\/\/6V\/ve975nL774Yg5rK94if\/nLX9x7RhVwfJruv\/9+u+SSS+zb3\/62a\/b+\/fvtW9\/6ljvg+zApoLVp0yZ5YNT\/10lKJ1tfpvvuu88UFtSDvGDBAl+a7QLBH\/\/4R3vooYe8aXNqQ9XzrR9WCpma1AuuXgYFfJ8m9Ryr1\/jee++1r371q140\/Wc\/+5n7IahjvaYdO3a4Ovrd737nRftvuOEGmzx5cvLHyNq1a23p0qX261\/\/Ovbt1xUq1X2it1K9x0OHDq30hIWRI0da3759beDAgbHfnkwbSEDLVCplvm7durlXTaj3LDHpbxs3brTmzZvnsMbiLKKANmPGDK96QCSlSzyf\/OQn3eVZTfplqHCsbfFt0olq6tSp7pKbgr8Pkw7oMldQVk\/gY4895kOzXRv1K1v1o8uCCva6vOzTI3QUDtT7p57L1157zR2D2rVr541\/oqGqIb33WL1Qvkz6IauhLfq+ajiCtuHVV19Nvrs57tuhkCPzc8891zX17bffdiHZhx+2ars6QhLDQHT80baogyEx6YeuvtO+9GhmUi8EtEyUUubRCVWP3VCXauqBXWleJyyfHsnha0BL7A7ti8cff9wWLlzofgmqV82X6aOPPrK7777b1VFZWZlrf2r3fVy3Q+3VL1UdID\/88EN3gPQpoOmyjkLaNddc4wK+xgLq2Ye+PIi6T58+Lpxp0klLPYLqSfbp8rgu62s7fvOb37i3tPgy6XK4LnH+7W9\/s1atWrkeZPWeffazn\/ViE1QjvXv3dvWisX96hJSOneoJjHvHQtWANm\/ePDeOTuMXE5OGvihE6zJuqUwEtCz3pC7paNzBK6+8UimgqftYJyuNC\/Fl8jmgqRdEQUE9CDrBJnrTfLFXuHzvvffsnXfecQ9G1gEo9WATx+2QuX6I\/OIXv3A\/RDQ+RDWvg6W+Fz4Mltb3tmnTpm6IgiYd1NXz6sNlHrVXA6F79uxpDzzwgGu\/xhLpBgdtR+LGhzjWTmqbNLB78+bNXvWeqf0agiBnjXlSDT377LOux0ZjMBs2bBh3dtdWjXfVf2vslsaQagyyLpHHvf1VA5p+0Cogpx4zdQVLvWg+9crWVjQEtNqE0vy7fomot+yyyy5z\/6pfVrqTcPv27V4FBV8D2okTJ0yv4NJ+8Onuu0Qp6eCYGigrKirsm9\/8pqufOB8o1fOkk5ROTprUE\/KPf\/zD9SDo17jupPJt0hguBR5duoqzfcI13TgbDa\/QvlFPrA+T7h7UdvjSaynTjz\/+2B3vNVb04osvdszqwVHN68HmHTp08IHetVG99xoicujQIdebppAW96lqQFNQ1jCj2bNnJ5uuY5OOrb4MFcnEnICWiVKVeTQoWndy6teIJl1m06Uf35K7bhmfOXOmd2PQ5P6vf\/0r2YugfaAue90dFvdHVaidvXr1su9\/\/\/sulKnnSY99UE35dpOJQo16n3y6SeDBBx90dwyOGjXKBTId4LUdvtyJqpPpr371K3cZR2OJVDMaaK\/xOD4ETJ1ANe6v6h15ORyG63wR3T2rmzHUi6zjjMZu6X+njuuq80Zl8YFTpkyxz3zmM+4yrY6VuhtVHQs+XB5XQFOPfeIO4KNHj7obAtR7r8ed6EeunqKgY5FPV7Fq230EtNqE0vy7fjnpgK6DpXoQdLlEA0c1LsGnSb9CdKD37S5O\/QLXXTyJnhztD3V362Tlw3gQHWh0iUq1o8tSehaRHtni0\/hF1blOsjpA+hTQ1GOm7+7zzz\/vakYH95\/85Cd23nnnefPV1Tg6jR3S8UbfAd2Zd8UVV3jRfo13krcuRfk2aZD6D3\/4w+T3Vj1pCgUa3uLDpOEU+kGiy\/yadPVBA+rjfmlcj\/Tp0aOHu2MzdaycOkX0PdDdnQr+GupSSndwah8R0CJ8s9TtrXDg011gETaXRfMsoAOPej3ifoDM82bHYnWJh6T68PytdGAaw6hjD6+jq\/ty0pAW1Y+v9qobjRdN9EbVvWB+P7GUX8tIQMtvrbA2BBBAAAEEEEAgsgABLTIhK0AAAQQQQAABBPIrQEDLrydrQwABBBBAAIE8CmjM4uuvv17pxrA8rj62qyKgxXbX0DAEEEAAAQT8F9C4Nz0zUWPfdFOd7sBMnXQn9bp169xd7QMGDHDPhUxMuqlKN7Lp0UQTJkzwHyOLLSCgZYHFrAgggAACCCCQnYBuqtCLzPWAYt1pn\/o6Jj2ySnex6z968LjeM6uHAesOWd1prZe76+5ZPV6GgJadO3MjgAACCCCAAAK1Cihk6Q5kvU82MenBsupVSzyPTY8z0TtnFebuvPNO9zDdFi1auAfT6vV4Pj0Sp1aQWmagBy2qIMsjgAACCCCAQK0C6QKa3oSh8KXes8Skv+mS529\/+1v3GCK9tUTPQbvtttvce1xDmQhooexpthMBBBBAAIECCPz973+3xYsXu7cTJN7J+9Zbb7kHEqe+eqlqQFNvmh7QrYfOpj5PVG9o0OsUEw\/v1gOmucRZgB3HKhFAAAEEEECgdAUUtDROTOHs4Ycfdi9k12VMvaP3yiuvTG541YCmmwK6dOni3m6QGtD0dgbdVFBKr23KZe\/Tg5aLGssggAACCCCAQFJANwJoQL\/ek6lHYihgaWxZ6pTuEqfGnqm3TC+j16Q3Negdodu3b3d3boY8EdBC3vtsOwIIIIAAAnkS0OMyBg0aZL169bI5c+a4V9mlTnpvr16RePvttyf\/rL\/pTs5p06a5vz3++OPukqeWD30ioIVeAWw\/AggggAACEQUUsvT4jBkzZtizzz5r7733nk2dOrXSu4bnz5\/vesjKy8uTn6ZnpD3xxBO2fv169yJ69bppuVatWkVskf+LE9D834dsAQIIIIAAAkUT0E0CI0aMsLlz51r79u1dOxTUNLh\/0qRJGbVLPWsKa6lj0TJasIRnIqCV8M5l0xBAAAEEECi0gAb7Hzp0yM4999xKH\/XBBx+4B9My5SZAQMvNjaUQQAABBBBAAIGCCRDQCkbLihFAAAEEEEAAgdwECGi5ubEUAggggAACCCBQMAECWsFoWTECCCCAAAIIIJCbAAEtNzeWQgABBBBAAAEECiZAQCsYLStGAAEEEEAAAQRyEyCg5ebGUggggAACCCCAQMEECGgFo2XFCCCAAAIIIIBAbgIEtNzcWAoBBBBAAAEEECiYAAGtYLSsGAEEEEAAAQQQyE2AgJabG0shgAACCCCAAAIFEyCgFYyWFSOAAAIIIIAAArkJENByc2MpBBBAAAEEEECgYAIEtILRsmIEEEAAAQQQQCA3AQJabm4shQACCCCAAAIIFEyAgFYwWlaMAAIIIIAAAgjkJkBAy82NpRBAAAEEEEAAgYIJENAKRsuKEUAAAQQQQACB3AQIaLm5sRQCCCCAAAIIIFAwgf8HpGp7\/mjQ0v8AAAAASUVORK5CYII=","height":419,"width":559}}
%---
%[output:14f0391e]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":2,"type":"double","value":[["0.480187691555994","0.480187691555994","0.480187691555994"],["0.480926504590384","0.480926504590384","0.480926504590384"]]}}
%---
%[output:4aa62d23]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":2,"type":"double","value":[["0.831709478944194"],["0.832989140657051"]]}}
%---
