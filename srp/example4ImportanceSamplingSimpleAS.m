%[text] # example for importance sampling of rendering equation (or SRP Ashikhmin–Shirley model)
%[text] モンテカルロ積分に使う乱数の数は `nMCint`
%[text] `それを nMC` 回行って平均と分散を計算し，plotする．（ややこしい）
%[text] 一応，SRPの計算にも対応（ただし太陽定数とかをかけていない）
clc
clear
cls

format long
%[text] ## Monte-Carlo condition & pre-allocation
nMC = 10;

srpTrue = zeros(nMC, 3);
srpImp = zeros(nMC, 3);
%[text] ## Ashikhmin–Shirley model
sat.F0 = 0.5 %[output:4ad78e7e]
sat.mCT = 0.3;

%[text] ### normal vector
sat.normal = [ 0 0 1 ];
% sat.normal = rand(1,3);
% sat.normal = sat.normal ./ vecnorm(sat.normal,2,2)
thetaN = deg2rad(40);
sat.normal = [0 sin(thetaN) cos(thetaN)];
sat.nu = 100;
sat.nv = 800;
sat.uu = [1 0 0];
sat.uv = cross(sat.normal, sat.uu);

dcm = triad(sat.uu, sat.normal, [1 0 0], [0 0 1]);
sat.qlb = dcm2q(4, dcm);

sat.rho = 0.5;
sat.area = 1;

thetaI = deg2rad(60);
phiI = deg2rad(20);

sunB = [sin(thetaI)*cos(phiI), sin(thetaI)*sin(phiI), cos(thetaI)];
sunB = sunB ./ norm(sunB);


% Monte-Carlo integration condition
nMCintUni = [100, 500, 1000, 5000, 10000, 50000, 100000];
nMCintUni = [10^5 10^6];

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

        NH = sat.normal * h'; % nFacet x nMC
        NS = sat.normal * sunB'; % nFacetx1
        VH = dot(v, h, 2)'; % 1xnMC
        NV = sat.normal * v'; % nFacet x nMC, 列方向に各reference vector, vの値

        F = sat.F0 + (1 - sat.F0) .* (1 - VH).^5; % nFacet x nMC
        k1 = sqrt((sat.nu+1) .* (sat.nv+1)) / 8 / pi; % nFacet x 1
        M = 1 ./ VH ./ max(NS, NV); % nFacet x nMC
        k2 = (sat.nu .* (sat.uu * h').^2 + sat.nv .* (sat.uv * h').^2) ./ (1 - NH.^2);
        D = NH .^ k2; % nFacet x nMC

        tmp = k1 .* F .* M .* D .* NV .* sin(thetaR)'; % nFacet x nMC

        csASx = NS .* (NS > 0) .* (NV > 0) .* tmp .* v(:,1)'; % nFacet x nMC matrix
        csASy = NS .* (NS > 0) .* (NV > 0) .* tmp .* v(:,2)';
        csASz = NS .* (NS > 0) .* (NV > 0) .* tmp .* v(:,3)';
      
        % % % for rendering eq.
        % csASx = (NS > 0) .* (NV > 0) .* tmp;
        % csASy = (NS > 0) .* (NV > 0) .* tmp;
        % csASz = (NS > 0) .* (NV > 0) .* tmp;

        % probability
        pdfAS = 1 / pi * 1 / (2 * pi); % scalar

        % specular component, nx3
        srpCs = [sum(csASx, 2), sum(csASy, 2), sum(csASz, 2)];
        srpCs = sum(srpCs,1) ./ pdfAS ./ N;

        srpTrue(j,:) = srpCs;
    end
    stdTrue(i,:) = std(srpTrue);
    meanTrue(i,:) = mean(srpTrue);
end
toc %[output:46d9a566]
%[text] ## 
%[text] ## importance sampling
% Monte-Carlo integration condition
% nMCint = 10:200:10000;
nMCintImp = [100 40^4];
nMCintImp = nMCintUni;

tic
for i = 1:length(nMCintImp)
    N = nMCintImp(i);
    for j = 1:nMC

%[text] #### transform to local frame
%[text] 法線ベクトルが${\\bf n} =\[0, 0, 1\]^T\n$になるlocal frameへ変換して計算する
        nFacet = size(sat.normal,1); % # of facet

        % sun vector and normal vectors at local frame (normal vector is along z-axis)
        sLocal = qRotation(4, repmat(sunB,nFacet,1), sat.qlb); % nFacet x 3
        nLocal = qRotation(4, sat.normal, sat.qlb); % (for debug), nFacet x 3

        NS = sat.normal * sunB'; % nFacetx1
%[text] ## sampling
%[text] importance sampling for ${\\bf{h}} = \[\\cos\\phi\_h \\sin\\theta\_h, \\sin\\phi\_h \\sin\\theta\_h, \\cos\\theta\_h\]^T$
        u1 = rand(N,1);
        u2 = rand(N,1);

        % phiH = atan(sqrt((sat.nu + 1) / (sat.nv + 1)) * tan(2*pi .* u1));
        % phiH = phiH .* (u1 <= 0.25) + (phiH + pi) .* (u1 > 0.25) .* (u1 < 0.75)...
        %     + (phiH + 2 * pi) .* (u1 >= 0.75);
        
        phiFun = @(u, nu, nv) atan(sqrt((nu + 1) / (nv + 1)) * tan(pi .* u / 2));
        ind1 = (u1 < 0.25);
        ind2 = (u1 >= 0.25) .* (u1 < 0.5);
        ind3 = (u1 >= 0.5) .* (u1 < 0.75);
        ind4 = (u1 >= 0.75);
        
        u1 = ind1 .* 4 .* u1 +  ind2 .* (4 * (u1 -0.25)) ...
            +  ind3 .* (4 * (u1 - 0.5)) + ind4 .* (4 * (u1 - 0.75));

        phiH = phiFun(u1, sat.nu, sat.nv); %一旦,計算
        phiH = ind1 .* phiH  + ind2 .* (phiH + pi/2) ...
            + ind3 .* (phiH + pi) + ind4 .* (phiH + 3/2*pi);            
        
        thetaH = acos((1-u2).^(1 ./ (sat.nu .* cos(phiH).^2 + sat.nv .* sin(phiH).^2 + 1)));

        h = [sin(thetaH).*cos(phiH) sin(thetaH).*sin(phiH) cos(thetaH)]; %Nx3

%[text] ### integration
        vx = 2 * (sLocal * h') .* h(:,1)' - sLocal(:,1); % nFacet x N
        vy = 2 * (sLocal * h') .* h(:,2)' - sLocal(:,2); % nFacet x N
        vz = 2 * (sLocal * h') .* h(:,3)' - sLocal(:,3); % nFacet x N

        VH = vx .* h(:,1)' + vy .* h(:,2)' + vz .* h(:,3)'; % nFacet x N
        NV = vz; % nFacet x N
        SH = sLocal * h'; % nFacet x N
        NH = h(:,3)'; % nFacet x N

        F = sat.F0 + (1 - sat.F0) .* (1 - VH).^5; % nFacet x nMC

        % weight
        W = abs(SH) .* NV ./ VH ./ NH./ max(NS, NV);        

        temp = W .* F; % nxN

        % for SRP
        csASx = NS .* (NS > 0) .*  (NV > 0) .* temp .* vx; % nFacetxN matrix
        csASy = NS .* (NS > 0) .*  (NV > 0) .* temp .* vy; % nFacetxN matrix
        csASz = NS .* (NS > 0) .*  (NV > 0) .* temp .* vz; % nFacetxN matrix

        % for rendering
        % csASx = (NS > 0) .* (NV > 0) .* temp; % nFacetxN matrix
        % csASy = (NS > 0) .* (NV > 0) .* temp; % nFacetxN matrix
        % csASz = (NS > 0) .* (NV > 0) .* temp; % nFacetxN matrix

        % specular component, nx3
        tmp = [mean(csASx,2), mean(csASy,2), mean(csASz,2)];
        tmp = qRotation(4, tmp, qInv(4, sat.qlb)); % transform back to the original body-fixed frame, nFacet x 3        

        sunlitFlag = (NS > 0); % nFacet x 1 matrix, 1: sunlit, 0: shade
        tmp = sat.area .*  sunlitFlag .* tmp;
        srpImp(j,:) = sum(tmp, 1);
    end

    stdImp(i,:) = std(srpImp);
    meanImp(i,:) = mean(srpImp);
end

toc %[output:844e8f0a]

%[text] ## show FIgs
figure %[output:23d7286a]
tiledlayout(3,1), nexttile %[output:23d7286a]
plotStd(nMCintUni, meanTrue(:,1), stdTrue(:,1), 'r'), hold on %[output:23d7286a]
plotStd(nMCintImp, meanImp(:,1), stdImp(:,1), 'b', 'b'); %[output:23d7286a]

nexttile %[output:23d7286a]
plotStd(nMCintUni, meanTrue(:,2), stdTrue(:,2)); hold on %[output:23d7286a]
plotStd(nMCintImp, meanImp(:,2), stdImp(:,2),'b', 'b'); %[output:23d7286a]

nexttile %[output:23d7286a]
plotStd(nMCintUni, meanTrue(:,3), stdTrue(:,3)); hold on %[output:23d7286a]
plotStd(nMCintImp, meanImp(:,3), stdImp(:,3), 'b', 'b'); %[output:23d7286a]

%[text] ### mean and standard deviation
[meanTrue(end,:) %[output:group:2343f67a] %[output:2148efe7]
    meanImp(end,:)] %[output:group:2343f67a] %[output:2148efe7]

[norm(meanTrue(end,:)) %[output:group:63f9bc25] %[output:22350528]
    norm(meanImp(end,:))] %[output:group:63f9bc25] %[output:22350528]
%[text] ## 
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40.7}
%---
%[output:4ad78e7e]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"sat","value":"    F0: 0.500000000000000\n"}}
%---
%[output:46d9a566]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 1.312812 秒です。\n","truncated":false}}
%---
%[output:844e8f0a]
%   data: {"dataType":"text","outputData":{"text":"経過時間は 1.364322 秒です。\n","truncated":false}}
%---
%[output:23d7286a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAmgAAAHOCAYAAAAsQvUnAAAAAXNSR0IArs4c6QAAIABJREFUeF7t3QtwVdW9+PEfoKQoj9LUVi0VKooKhYHCFShYHsoA9Tm2MCmDDThoKRAdoOVRitfWBkVhYiRJeUgFopQSEFsLFJACEqFkKqGICKIZag3KozhIkT\/v\/\/zWvefcQ0jCzlnnsfde3z3TaRv22mevz\/qdvX5n7bXXrnfhwoULwoYAAggggAACCCDgG4F6JGi+aQtOBAEEEEAAAQQQMAIkaAQCAggggAACCCDgMwESNJ81CKeDAAIIIIAAAgiQoBEDCCCAAAIIIICAzwRI0HzWIJwOAggggAACCCBAgkYMIIAAAggggAACPhMgQfNZg3A6CCCAAAIIIIAACRoxgAACCCCAAAII+EyABM1nDcLpIIAAAggggAACJGjEAAIIIIAAAggg4DMBEjSfNQingwACCCCAAAIIkKARAwgggAACCCCAgM8ESNB81iCcDgIIIIAAAgggQIJGDCCAAAIIIIAAAj4TIEHzWYNwOggggAACCCCAAAkaMYAAAggggAACCPhMgATNZw3C6SCAAAIIIIAAAiRoxAACCCCAAAIIIOAzARI0nzUIp4MAAggggAACCJCgEQMIIIAAAggggIDPBEjQfNYgnA4CCCCAAAIIIECCRgwggAACCCCAAAI+EyBB81mDcDoIIIAAAggggAAJGjGAAAIIIIAAAgj4TIAEzWcNwukggAACCCCAAAIkaMQAAggggAACCCDgMwESNJ81CKeDAAIIIIAAAgiQoBEDCCCAAAIIIICAzwRI0HzWIJwOAggggAACCCBAgkYMIIAAAggggAACPhMgQfNZg3A6CCCAAAIIIIAACVqVGDh48KCsWLFCKisrZcCAAdKtWzdp0KABkYIAAggggAACCKRMgAQthvrYsWMyePBg6dmzp\/Tt21eWLl0q9evXl7y8vJQ1CB+EAAIIIIAAAgiQoMXEQElJiZSWlkp+fr7566lTp6RXr16ybNkyadGiRbXR8tBDD0lZWRmRhAACCCCAAAJJFrj99tuluLg4yZ\/ij8OToMW0w6RJk6Rjx46SlZUV\/evkyZPNrU5N1KrbbrnlFtm7d68\/WjNBZ0GdEgSZ5MPQTkkGTtDhaacEQSb5MLRTkoETdPgwtlNNNCRoMTLZ2dkycuRI6d69e\/Svs2bNkqZNm4r+Gwlagr5haThMGL\/U1CkNgRTHR9JOcaCloQjtlAb0OD4yjO1EguYhEB5\/\/HEzehaboOXm5krLli1l6NChJGgeDP26Sxi\/1NTJr9F28XnRTrRTugSIvXTJJ+ZzGUGLcSwsLJRmzZpdlIwNGzbMjJ716dOHBC0xMZeWo3ChSgt7nT+UdqozWVoK0E5pYa\/zh9JOdSbzVQEStJjmKC8vF52HtnjxYsnMzJQtW7bIuHHjZNOmTZKRkeFMglZQUCBjxozxVaDangx1shVMTXnaKTXOtp9CO9kKpqZ8GNspjElnTdFAglZFZvXq1bJkyRKpqKiQxo0by4wZM6Rdu3Y1fptcCpbUXFL4FAQQQAABBKoXcKnPJUGr4Vtw8uRJadSo0WW\/Iy4Fy2Ux2AEBBBBAAIEkCrjU55KgWQaSS8FiSUVxBBBAAAEErARc6nNJ0KxCRcSlYLGkojgCCCCAAAJWAi71uSRoVqFCgmbJR3EEEEAAAQQ8C5CgeaZiR5eChdZGAAEEEEAgnQIu9bmMoFlGmkvBYklFcQQQQAABBKwEXOpzSdCsQoVbnJZ8FEcAAQQQQMCzAAmaZyp2dClYaG0EEEAAAQTSKeBSn8sImmWkuRQsllQURwABBBBAwErApT6XBM0qVLjFaclHcQQQQAABBDwLkKB5pmJHl4KF1kYAAQQQQCCdAi71uYygWUaaS8FiSUVxBBBAAAEErARc6nNJ0KxChVuclnwURwABBBBAwLMACZpnKn\/tePDgQVmxYoVUVlbKgAEDpFu3btKgQYNLTvLdd9+VVatWyYULF+T+++83r2uKbDNnzoyWOXfunHlh+qhRo2qsqEvB4q\/W5mwQQAABBFwTcKnPDc0I2rFjx2Tw4MHSs2dP6du3ryxdulTq168veXl5F8Xvrl27ZMSIEeY\/N998szzzzDOSnZ0tWVlZoglZ27ZtpaioyJQ5f\/681KtXT+666y4SNNeuAtQXAQQQQMB3AiRovmuSy59QSUmJlJaWSn5+vtn51KlT0qtXL1m2bJm0aNEieoAJEyZIhw4dZOjQoeZve\/bskUceeUQ2b94sR48elUGDBsn69esv\/4H\/u0fs6Ft1hcaMGSM5OTmej8eOCCCAAAIIuCpQUFAgs2bNqrX6e\/fudYInNCNokyZNko4dO5qRsMg2efJkc6tTE7XIprc9i4uLzehZ7N\/WrVsneot06tSpMnfuXNm+fbs5XrNmzWoNBJeyeSe+EVQSAQQQQMC3Ai71uaFJ0PQ25ciRI6V79+7RwNIsvGnTpuYWpm5nz56Vdu3ayc6dOyUjIyO635AhQ2TKlCmit0mHDx8uAwcOlNatW4smbf369at1BMylYPHtN5YTQwABBBBwQsClPjeQCZpO7h8\/frw0bNjQzBO7++675dVXXzWjZ7EJWm5urrRs2TJ6O1PLde7cWbZu3XpRgqYJmc5V0+Pt27dP+vfvbwL9+PHjZvRt48aNJtGrbnMpWJz49lNJBBBAAAHfCrjU5wYyQdPIKSsrkyuvvNIkaDrHTOea6e3IyNwy3WfYsGFm9KxPnz7RYNN\/19Gy2267zfzt9OnT0r59e3M8La9JnD4YENk0WZs2bZpJ7EjQfPud5cQQQAABBBwQIEELYCOXl5eLzkNbvHixZGZmypYtW2TcuHGyadMm88DA\/v37zcMB8+bNE32SU5fT0G3+\/PnmlmdhYaGsXLnSPP2pCdl1111nHhaYOHGiOUaTJk1I0AIYF5wyAggggEB4BEjQAtqWq1evliVLlkhFRYU0btxYZsyYYeacbdu2TZ544glZs2aNnDlzRhYsWGD+98cff2yStunTp0vz5s3NHDV9gEDXUjt06JAZUdPbpF26dKlRxKVgCWhYcNoIIIAAAiERcKnPDewtztpi7eTJk2aB2do2vTWqyVrswwKx+3s5hu7vUrCE5PtNNRBAAAEEAirgUp8bygQtlXHnUrCk0pXPQgABBBBAoKqAS30uCZpl\/LsULJZUFEcAAQQQQMBKwKU+lwTNKlS4xWnJR3EEEEAAAQQ8C5CgeaZiR5eChdZGAAEEEEAgnQIu9bmMoFlGmkvBYklFcQQQQAABBKwEXOpzSdCsQoVbnJZ8FEcAAQQQQMCzAAmaZyp2dClYaG0EEEAAAQTSKeBSn8sImmWkuRQsllQURwABBBBAwErApT6XBM0qVLjFaclHcQQQQAABBDwLkKB5pmJHl4KF1kYAAQQQQCCdAi71uYygWUaaS8FiSUVxBBBAAAEErARc6nNJ0KxChVuclnwURwABBBBAwLMACZpnquDsePr0aSkpKZFdu3ZJ165dpX\/\/\/tW+UH3Dhg2yfft2OXfunOTk5Fz2pesuBUtwWpszRQABBBAIo4BLfa4zI2hjx46VY8eOSXZ2tuzcuVOWLl0qa9euvSQBe++99+TAgQMyatQoefvtt6Vx48a1xrhLwRLGLzt1QgABBBAIjoBLfa4TCVpFRYUMGTJENm3aJBkZGSYSR48eLf369ZMHHnig2sjUINi9e7c0aNDgsglabTuMGTPGjMSxIYAAAggggEDtAgUFBTJr1qxad9q7d68TjE4kaMuXL5fS0lLJy8uLNuprr70mH374oYwfP77GBG3Pnj1Sr169yyZorgSLE98IKokAAggg4FsBRtB82zTxndjs2bPlzJkzF41kbdu2TRYuXChFRUUkaPGxUgoBBBBAAIGUCpCgpZQ78R9WWFgo\/\/rXv8yB27RpIw0bNpTPPvvsogRt3bp1oqNoum91mwYBI2iJbxuOiAACCCCAQLwCJGjxyvmknCZWJ0+eNGejk\/yPHj0qxcXFove2I9ucOXPMQwMTJkwgQfNJu3EaCCCAAAII1CZAghay+Dh+\/Lh5IEBvdXbs2FEqKyslKytL5s2bZ0bYysvLpW3bthc90alBoE901q9fv1YNl4IlZGFBdRBAAAEEAibgUp\/rxEMCGn+6tMZLL71knszUkbNJkyaZJzh1pK1Hjx7y6quvSqtWrUyoRv6mT302adKEBC1gX2BOFwEEEEAgnAIkaOFs12jy1ahRo4TV0KVgSRgaB0IAAQQQQCAOAZf6XGdG0OKIA09FXAoWTyDshAACCCCAQJIEXOpzSdAsg8ilYLGkojgCCCCAAAJWAi71uSRoVqHCy9It+SiOAAIIIICAZwESNM9U7OhSsNDaCCCAAAIIpFPApT6XETTLSHMpWCypKI4AAggggICVgEt9LgmaVahwi9OSj+IIIIAAAgh4FiBB80zFji4FC62NAAIIIIBAOgVc6nMZQbOMNJeCxZKK4ggggAACCFgJuNTnkqBZhQq3OC35KI4AAggggIBnARI0z1Ts6FKw0NoIIIAAAgikU8ClPpcRNMtIcylYLKkojgACCCCAgJWAS30uCZpVqHCL05KP4ggggAACCHgWIEHzTOWvHQ8ePCgrVqyQyspKGTBggHTr1k0aNGhwyUkuWrRIDh8+LM2aNZMRI0Zc9O8zZ86Mljl37pzoi9VHjRpVY0VdChZ\/tTZngwACCCDgmoBLfW5oRtCOHTsmgwcPlp49e0rfvn1l6dKlUr9+fcnLy7skft966y3ZvXu3rFmzRpYtWxb9d03I2rZtK0VFReZv58+fl3r16sldd91FgubaVYD6IoAAAgj4ToAEzXdNcvkTKikpkdLSUsnPzzc7nzp1Snr16mUSsBYtWlxygB07dpjkbeHChdF\/O3r0qAwaNEjWr19\/+Q\/83z00WGrbxowZIzk5OZ6Px44IIIAAAgi4KlBQUCCzZs2qtfp79+51gic0I2iTJk2Sjh07SlZWVrThJk+ebG51aqJWddME7fnnn5cFCxZE\/+mDDz6QqVOnyty5c2X79u3meHobtLbNpWzeiW8ElUQAAQQQ8K2AS31uaBK07OxsGTlypHTv3j0aWJqFN23aVPTfvCRoW7ZskeHDh8vAgQOldevWsm7dOunXr1+tI2AuBYtvv7GcGAIIIICAEwIu9bmBTNAuXLgg48ePl4YNG5p5Ynfffbe8+uqrZvQsNkHLzc2Vli1bytChQz0laBUVFbJv3z7p37+\/2f\/48eNm9G3jxo0m0atucylYnPj2U0kEEEAAAd8KuNTnBjJB08gpKyuTK6+80iRoOsdM55rp7cjYZGzYsGFm9KxPnz6eEjTdSZM\/fTAgsmmyNm3aNOncuTMJmm+\/spwYAggggIALAiRoAWzl8vJy0XloixcvlszMTNHblePGjZNNmzaZBwb2798vHTp0iNZM99cHCmLnoK1cudI8\/akJ2XXXXWceFpg4caI5RpMmTUjQAhgXnDICCCCAQHgESNAC2parV6+WJUuWiN6qbNy4scyYMUPatWsn27ZtkyeeeMIsqxHZ9G+6nEbsU5xnz56V4uJis5baoUOHzIic3ibt0qVLjSIuBUtAw4LTRgABBBAIiYBLfW5gb3HWFmsnT540C8zabF6P4VKw2HhSFgEEEEAAAVsBl\/rcUCZotgFQl\/IuBUtdXNgXAQQQQACBRAu41OeSoFlGj0vBYklFcQQQQAABBKwEXOpzSdCsQoWXpVvyURwBBBBAAAHPAiRonqnY0aVgobURQAABBBBIp4BLfS4jaJaR5lKwWFJRHAEEEEAAASsBl\/pcEjSrUOEWpyUfxRFAAAEEEPAsQILmmYodXQoWWhsBBBBAAIF0CrjU5zKCZhlpLgWLJRXFEUAAAQQQsBJwqc8lQbMKFW5xWvJRHAEEEEAAAc8CJGieqdjRpWChtRFAAAEEEEingEt9LiNolpHmUrBYUlEcAQQQQAABKwGX+lwSNKtQEfnqV78q7du3tzzK\/xXv3bt3wo7lyoF69erlSlUTWk9iLaGcHAwBBFIgQIKWAuRUf8Tp06elpKREdu3aJV27dpX+\/ftX+0L1AwcOyCuvvCKffPKJtGnTRu677z65\/vrrazzdFi1aSH5+fsKqs3Xr1oQdy5UDbdmyxZWqJrSexFrdOUlq626mJXCruxs\/PKs3+\/GPfywfffRR3UEDWMKZEbSxY8fKsWPHJDs7W3bu3ClLly6VtWvXXpSkffbZZ3LnnXfKo48+Kv\/1X\/8l69atk5UrV8qGDRvkiiuuqLZ5b7rpJnnzzTcD2PScMgII1FWApLauYv+zP251d3Pqh2cdBiauuuoqeePEibqDBrCEEwlaRUWFDBkyRDZt2iQZGRmmmUaPHi39+vWTBx54INpsr732mknaioqKzN\/OnTsnPXr0kEWLFpnRtOq2QS1ayIULF2ps+uHDh0unTp0SExp1COLEfCBHCa0Ao46hbdqUVoxrUkq5Q\/th3btfVLXyHTsuqer\/O3nS\/G3HVVfJUyRo4QmF5cuXS2lpqeTl5V2UjH344Ycyfvz46N8+\/\/xzOXHihFx33XXmb5HEbvPmzXLllVdWC5KbmSnDb7klNVjf\/W5KPud3v\/udPPzwwyn5rFR9CHWqIl3lgpiqdrjc59BOlxPyx79H26lKHPl1pMyv5+XHUTK\/WkVuk7\/zzjty5MgRf3wRknwWToygzZ49W86cOSM5OTlRzm3btsnChQujo2VVndevXy9TpkyR\/\/7v\/5aBAwfW2AxhvMV5xx13iCalYdqo08Wt6ceLsJ6THxM02050x44d0rFjx4R+nfzYflpBv84183Jes2bNuqiPSGiD1XCwZM8z0\/laegeoLpsXq7ocL9H78pBAokVTfLzCwkL517\/+ZT5Vb002bNhQdH5ZbIKm88v0lqbuG7tpIjd16lTR0bXp06fLjTfeWOvZJzNBS9dF+LHHHpMXXnih2nqn65y8hFBtHWkyOkkv5+RXLz9ehPWc0tFJXq4dbTvReDrJy51TutsvjJ0kdbpc1Pnj38PYTjXJhnIEbc+ePXLyf+9XN27cWI4ePSrFxcVSUFAQdZgzZ455aGDChAkX2Tz99NPm\/\/\/85z+v8cGA2AJXX321fPHFF0mJ3HRdhMvKyuT222+vtk7pOicvwLV1pMnoJL2cUzK9wnihok5eoir9+9BO6W8DL2dAO3lR8u8+oUzQqnIfP37cPBCgtzr1VkNlZaVkZWXJvHnzzAhbeXm5tG3bVjQxeeaZZ+TPf\/6zNGjQwBxGHwA4f\/589P9XPbaugxY7MuffpvZ+Zn4cxfB+9tXvSZ1sBVNTnnZKjbPtp9BOtoKpKR\/Gdpo\/fz7LbKQmfFL3Kbq0xksvvSS7d+82I2eTJk0yT3DqSJs+qfnqq6+apTeWLFkiOiqmm97u1Fujzz\/\/fI3z0L7\/\/e+bp0PZEEAAAQQQQCC5AnqnZNWqVcn9EJ8c3YkRtFhrTcgaNWrkE35OAwEEEEAAAQQQuFTAuQSNIEAAAQQQQAABBPwuQILm9xbi\/BBAAAEEEEDAOQESNOeanAojgAACCCCAgN8FSND83kKcHwIIIIAAAgg4J0CC5lyTU2EEEEAAAQQQ8LsACZrfW4jzQwABBBBAAAHnBEjQ4mxyfb\/Z4cOHpVmzZjJixIg4j+KvYm+\/\/bZZpFeXImnfvr08+OCDgV+SRN+p+te\/\/lXOnj1r3o6gdapXr56\/4OM8m\/fff1+0zX70ox\/FeQR\/FDtw4IC88sorcuWVV5oTOnfunNx6661y9913++MELc7i4MGDZn3FiooK6dy5s\/zgBz8I9HdKFz49ffp0dOFubSv9jy7WHeTli3Sxcr1OHDp0yCxmPnjw4BoXJ7cIh5QW3bBhg6mTruf5wx\/+ULp06ZLSz0\/Eh2kdtm\/fXm2MaRyWlJTIrl27pGvXrtK\/f\/9Ax2B1XiRocUbRW2+9ZRa9XbNmjSxbtizOo\/inmL48\/qc\/\/ank5uZK06ZNzYvkv\/zlL8uzzz7rn5Os45nou1Z\/+9vfmjqdOHFCZs6caS68Q4cOreOR\/Le7JtH33HOPfPOb35QFCxb47wTrcEaaZOobPEaOHGlK6YX32muvlU6dOtXhKP7bVRe5vu+++0y8acKpr5f71re+ZeIxqNvq1avNqdevX9\/8t3aeb7zxhujf9Z3HQdz27dsnDz\/8sPziF78w1z59P7Mm0+PHjw9idcw5L1++XIqKiuTJJ5+UTz75xLzjVq\/l3bt3D1Sd3nvvPdEfcKNGjTI\/RvXVjZFt7NixZtH57Oxs0YXo9YfQ2rVrQ5WkkaBZhKu+gDsvL88kM0HfpkyZIvri9+HDh5uq6MviBw0aZC7AQd00Qfv6178evSjp\/9eORDvKoG8TJ040b7nQV5HpK8uCvGkH\/+abb8qvf\/3rIFfjknPXUXb9EafJp2464q4jAvojIQybjkrriLS+t\/iOO+4IbJWee+4582P0kUceMXX4+9\/\/bmLxT3\/6U2DrNHDgQJk2bVr0R87KlSvl5Zdflt\/\/\/veBrJO+U1S\/S5FXMOqI9JAhQ8xbfDIyMkydRo8ebV7pqG8ICstGgmbRkpqg6Wuggj6CoQR6K+ZLX\/qSuWWrm\/4C27x5s6lfGDbtTKZPn25uo02YMCHQVdKLrLaPJtU6Qjh37txA10d\/+Wr86e0\/\/WGgt6IjF90gV0w7fB0V1FHOf\/zjH3LzzTdLq1atglyli85d47C0tNSMOAV50x\/Y7777rrk+6PQHrdff\/vY3yc\/PD2y1NKHRtrnmmmtMHfbv32+S6aD+4Nb67NmzJzo9Ra9\/Wj8dIIls+gNcrx9BHvmsGnAkaBZfwTAlaBEGTWT0ZbT63lL9xaWjakHePv\/8c3nsscfMEPgNN9xg6hQ7TB60umk99JeiXqD+\/e9\/mwtU0BM0vf2iSdr3vvc98wNB5w3qu3L79OkTtOa56HzvvPNOk5zpph2MjhTqCHUYbrHrlAGt3x\/+8Adp2bJloNtJb6nrLc4PPvhAmjdvbkamdfTsa1\/7WmDrpTHWt29fE286R3DGjBnmmq6jg02aNAlcvaomaLNnzzZz63TuY2TTaTqabOut3bBsJGgWLRm2BE1HMbTz11\/62kFGRtMsiNJeVBPOTz\/9VD766CNZvHix6Shjv9RpP8E6nIC2jw7rv\/DCC9KuXTvR+RmaoOnFSm91Rob\/63BIX+y6detWufrqq6VDhw7mfPRCqyO3Qb0dE0HVScs9e\/aUqVOnmj\/pXCB98EHrF3kgwhcNEMdJ6OTsjRs3Bn70TKuuUx60TXQuk8bh66+\/bkZi9M7IFVdcEYdO+ovo+eucW\/1vnael81V1vrTeYg9inaomaPpDWxPp2Gv5unXrREfRgj6iGxs9JGgW36UwJWinTp2SYcOGmV\/3YXh6TptVL0yxSWZlZaXce++9UlZWFsiLlI4yaWeinYhuOorx8ccfm1\/6+gtZn2QKw6ZztTSx0dtOQexMIm1Q3ZyYbt26mdFCHc0N8qZPBWr9gj7Kef78ebntttvM3NQbb7zRNImOzOh3acmSJdKmTZsgN5PoHQSdunLkyBEzmqZJWhC3qgmaJtTFxcVSUFAQrY5eG\/WaH\/QpLCRoCYpQfTRb5ymEYQ6a\/tr64osvor\/2lUiHxvVprSAuS6Hn3rt3b\/nZz35mkjIdYdKlHHRCvc6tC8OmCYyONAX9IYFf\/epX5gnAMWPGmIRML7pat6B\/r7QzfPHFF80tF50LpHGnE+p17kyQE0\/tBHWeYNWn6oL6ndInbfXBDR2d1mudztPS\/x07hytodXv66aflK1\/5irl1q9dwfUJVl04K6u11TdD0jkHk6eHjx4+bBwL07oEui6I\/vrOyssy1UJ+YDsvGCJpFS2oWrxffMDzFqb+I9cmYyOiM\/orUIWTtVII6F0O\/0Hp7SUeZ9JaSrgOka9bp7cEwbNpB6gUq6AmajphpMvaXv\/zFxJxecH\/zm9\/I9ddfH\/hm0vl1OvdH5zbpd0ufrPv2t78d6HrpPCZtH72dFIZNJ5\/\/8pe\/jF4ndCRNO3t9EjKom07p0B86On1AN70ropPng3hrXZcU6tGjh3liM3b+nM7H1e+WPt2pPxp0Wk6YnuDUdiNBC+o3kPP2LKBfcB2xCOLFyXMlQ7BjZNHToK6nVVMT6DxI\/cET5IVcQxBel62CPiygMRimdtK407mpkZGnyyIEcAe9voepzbjFGcAg5JQRQAABBBBwSUCnpkSm2OgPnSBPDYin3RhBi0eNMggggAACCCCQVAFdYPePf\/yjfOMb3zBzpHWKQFCfVo8HigQtHjXKIIAAAggggIAnAb3VqksCaXKly+noBP\/YTR8KWrVqlXmY6\/777zfLIemmy2no\/pEleDx9WIh2IkELUWNSFQQQQAABBPwmoHP79D2ZunaePnQWu9q\/vuxcH97S\/+ganPpqNF2TTh\/U0KfU9cEhLa8LIkdex+W3+iXrfEjQkiXLcRFAAAEEEEAgKqDLzuhcMn0FWmTTdct0hCyyBIg+VauJmK4goGtW6miarmepL3vX1QYi69W5wEqC5kIrU0cEEEAAAQTSLFBdgqaLN+uiszp6Ftn0b\/pmAB1d+853vmPezasvtR80aFCo3md7ueYgQbucEP+OAAIIIIAAAjUK\/POf\/5RFixaZBXEjk\/jff\/99s1Ze7Mr+VRM0HU3TdSl1TTNNwiKbLhQ8ZcoU+c9\/\/mNe+aZrCF577bWBfU1fvKFDghavHOUQQAABBBBAwNy21MV+NTl76qmnzDtA9TamvoKuU6dOUaGqCZo+FNC5c2ezoG5sgqaLBOtDBZG3Ari4xIaikaDx5UIAAQQQQAABKwGdyP\/kk0+KvobpnXfeMQlW1acvq7vFqXPPdLRM34mqmy4YrK+l0vlnse9Stjq5gBYmQQtow3HaCCCAAAII+ElAl8t48MEHzXuQCwsLL1lYVl9Lpy+o\/8lPfhI9bf2bzjXT90HrNn\/+fHPLU8u7vpGguR4B1B\/zUaGYAAAgAElEQVQBBBBAAAFLAU2ydPkMXRrj9ddfl08\/\/VSmT59+0Sv25syZY0bIcnJyop+ma6Tpu3jXrFlj3oeqo25aTt9f6\/pGguZ6BFB\/BBBAAAEELAT0IYFHH31UioqKpHXr1uZImqgdPnxYcnNzPR1ZR9Y0WYudi+apYIh3IkELceNSNQQQQAABBJItoJP9jxw5Itdcc81FH3Xo0CGzMC1bfAIkaPG5UQoBBBBAAAEEEEiaAAla0mg5MAIIIIAAAgggEJ8ACVp8bpRCAAEEEEAAAQSSJkCCljRaDowAAggggAACCMQnQIIWnxulEEAAAQQQQACBpAmQoCWNlgMjgAACCCCAAALxCZCgxedGKQQQQAABBBBAIGkCJGhJo+XACCCAAAIIIIBAfAIkaPG5UQoBBBBAAAEEEEiaAAla0mg5MAIIIIAAAgggEJ8ACVp8bpRCAAEEEEAAAQSSJkCCljRaDowAAggggAACCMQnQIIWnxulEEAAAQQQQACBpAmQoCWNlgMjgAACCCCAAALxCYQqQTt48KCsWLFCKisrZcCAAdKtWzdp0KDBJTKLFi2Sw4cPS7NmzWTEiBEX\/fvMmTOjZc6dOyeNGjWSUaNGxadLKQQQQAABBBBAIA6B0CRox44dk8GDB0vPnj2lb9++snTpUqlfv77k5eVdwvLWW2\/J7t27Zc2aNbJs2bLov2tC1rZtWykqKjJ\/O3\/+vNSrV0\/uuuuuOGgpggACCCCAAAIIxCcQmgStpKRESktLJT8\/30icOnVKevXqZRKwFi1aXKKzY8cOk7wtXLgw+m9Hjx6VQYMGyfr16+PTpBQCCCCAAAIIIJAAgdAkaJMmTZKOHTtKVlZWlGXy5MnmVqcmalU3TdCef\/55WbBgQfSfPvjgA5k6darMnTtXtm\/fbo6nt0Fr2x566CEpKytLQFNwCAQQQAABBBCoTeD222+X4uJiJ5BCk6BlZ2fLyJEjpXv37tGGmzVrljRt2lT037wkaFu2bJHhw4fLwIEDpXXr1rJu3Trp16+f5OTk1BgMt9xyi+zduzdUwUKdgtGctBPtlC4BYi9d8nX7XNqpbl5+2zs0Cdrjjz9uRs9iE7Tc3Fxp2bKlDB061FOCVlFRIfv27ZP+\/fub\/Y8fP25G3zZu3GgSveo2vgB+C+nqz4d2op3SJUDspUu+bp9LO9XNK117h7GdarIMTYJWWFhobkfGJmPDhg0zo2d9+vTxlKDpThcuXDAPBkQ2TdamTZsmnTt3JkFL1zcyAZ8bxi81dUpAYKTgELRTCpAT8BG0UwIQU3CIMLZT6BO08vJy0XloixcvlszMTNHblePGjZNNmzaZBwb2798vHTp0iDro\/vpAQewctJUrV5qnPzUhu+6668zDAhMnTjTHaNKkiTMJWkFBgYwZMyYFX7XUfQR1Sp21zSfRTjZ6qStLO6XO2uaTwthOJGg2EZHGsqtXr5YlS5aI3qps3LixzJgxQ9q1ayfbtm2TJ554wiyrEdn0b7qcRuxTnGfPnjWTD3UttUOHDpkROb1N2qVLlxpr5VKwpLFp+WgEEEAAAQTEpT43NLc4Y+P25MmTZoFZm83rMVwKFhtPyiKAAAIIIGAr4FKfG8oEzTYA6lLepWCpiwv7IoAAAgggkGgBl\/pcEjTL6HEpWCypKI4AAggggICVgEt9LgmaVaiIU\/fDLakojgACCCCAgJUACZoVn1uFXQoWt1qW2iKAAAII+E3ApT6XETTL6HMpWCypKI4AAggggICVgEt9LgmaVahwi9OSj+IIIIAAAgh4FiBB80zFji4FC62NAAIIIIBAOgVc6nMZQbOMNJeCxZKK4ggggAACCFgJuNTnkqBZhQq3OC35KI4AAggggIBnARI0z1Ts6FKw0NoIIIAAAgikU8ClPpcRNMtI++pXvyo5OTnRo\/Tq1cvyiMkv3rt37+R\/CJ+AAAIIIIBAggVI0BIMGubDXX311fLFF1+EuYopr1sQEsggnGPKG44PRAABBJIsMH\/+fPnoo4+S\/Cn+ODwjaJbt0KJFC8nKyrI8CsURSKzAli1bEnvAJBxt69atSTgqh\/SbQBB+zAThHP3Wruk6n1mzZsmRI0fS9fEp\/VwSNEvu3MxMGX7LLf9zlO9+1\/JoFA+cQPfugTtlTti7QI1JpI\/a3e+JLj8WvMdbTXsGIYFM1TkygmYfT84cYVCLFpKfn\/8\/9WVEwJl2j1Y0ACNV7jVKCmrMdz0FyD77iDQm5X5Pwk33l6Lm2nHVVfLGiRMp+rT0fgwjaJb+N910k7z55puWR6F4EASCcJEMgmO6z5F2THcLJO7zUzo6R1KeuIaLOVJdE7urrrpKTpCgJaUtQnfQzMxM0adKatvoEMLR7Kkawg+Hln9rQTv6t23qemZBeGq+rnVybf+6fh95itO1CLGo7w033CCLFi2q9Qh1DUCL06EoAggggAACoRUgQQtt0ya+Yi4FS+L1OCICCCCAAALeBVzqc5mD5j0uqt3TpWCxpKI4AggggAACVgIu9bkkaFahwrs4LfkojgACCCCAgGcBEjTPVMHZ8fTp01JSUiK7du2Srl27Sv\/+\/aVRo0aXVGDDhg2yfft2OXfunHmFU3X7xBZyKViC09qcKQIIIIBAGAVc6nOdGUEbO3asHDt2TLKzs2Xnzp2ydOlSWbt27SUJ2HvvvScHDhyQUaNGydtvvy2NGzeuNcZdCpYwftmpEwIIIIBAcARc6nOdSNAqKipkyJAhsmnTJsnIyDCROHr0aOnXr5888MAD1UamBsHu3bulQYMGJGjB+e5ypggggAACIRYgQQtZ4y5fvlxKS0slLy8vWrPXXntNPvzwQxk\/fnyNCdqePXukXr16l03QatthzJgx5lYpGwIIIIAAAgjULlBQUCD6vs3atr179zrB6MQI2uzZs+XMmTMXJUrbtm2ThQsXSlFRkXWC5kqwOPGNoJIIIIAAAr4VYATNt00T34m9\/PLL8tlnn12UoK1bt050FK2wsJAELT5WSiGAAAIIIJBSARK0lHIn\/8N0tKy4uFh06DSyzZkzxzw0MGHCBBK05DcBn4AAAggggIC1AAmaNaG\/DnD8+HHzQIDe6uzYsaNUVlZKVlaWzJs3T9q0aSPl5eXStm3bi57o1CDQJzrr169fa2VcChZ\/tSpngwACCCDgmoBLfa4Tc9A0gHVpjZdeesk8makjZ5MmTTJPcJ48eVJ69Oghr776qrRq1crEeuRv+tRnkyZNSNBcuwJQXwQQQAABXwqQoPmyWRJzUpp8XW7x2bp8kkvBUhcX9kUAAQQQQCDRAi71uc6MoCU6SCLHcylYkmXIcRFAAAEEEPAi4FKfS4LmJSJq2celYLGkojgCCCCAAAJWAi71uSRoVqHCy9It+SiOAAIIIICAZwESNM9U7OhSsNDaCCCAAAIIpFPApT6XETTLSHMpWCypKI4AAggggICVgEt9LgmaVahwi9OSj+IIIIAAAgh4FiBB80zFji4FC62NAAIIIIBAOgVc6nMZQbOMNJeCxZKK4ggggAACCFgJuNTnkqBZhQq3OC35KI4AAggggIBnARI0z1Ts6FKw0NoIIIAAAgikU8ClPpcRNMtIcylYLKkojgACCCCAgJWAS30uCZpVqHCL05KP4ggggAACCHgWIEHzTMWOLgULrY0AAggggEA6BVzqcxlBs4w0l4LFkoriCCCAAAIIWAm41OeSoFmFCrc4LfkojgACCCCAgGcBEjTPVOzoUrDQ2ggggAACCKRTwKU+lxE0y0hzKVgsqSiOAAIIIICAlYBLfS4JmlWocIvTko\/iCCCAAAIIeBYgQfNMxY4uBQutjQACCCCAQDoFXOpzGUGzjDSXgsWSiuIIIIAAAghYCbjU55KgWYUKtzgt+SiOAAIIIICAZwESNM9U\/trx4MGDsmLFCqmsrJQBAwZIt27dpEGDBpec5LvvviurVq2SCxcuyP333y\/a4JFt5syZ0TLnzp2TRo0ayahRo2qsqEvB4q\/W5mwQQAABBFwTcKnPDc0I2rFjx2Tw4MHSs2dP6du3ryxdulTq168veXl5F8Xvrl27ZMSIEeY\/N998szzzzDOSnZ0tWVlZoglZ27ZtpaioyJQ5f\/681KtXT+666y4SNNeuAtQXAQQQQMB3AiRovmuSy59QSUmJlJaWSn5+vtn51KlT0qtXL1m2bJm0aNEieoAJEyZIhw4dZOjQoeZve\/bskUceeUQ2b94sR48elUGDBsn69esv\/4H\/u4dLweIZhR0RQAABBBBIgoBLfW5oRtAmTZokHTt2NCNhkW3y5MnmVqcmapFNb3sWFxeb0bPYv61bt070FunUqVNl7ty5sn37dnO8Zs2a1RpisbdHq9txzJgxkpOTk4Qw5ZAIIIAAAgiES6CgoEBmzZpVa6X27t0brkrXUJvQJGh6m3LkyJHSvXv3aFW1kZs2bWpuYep29uxZadeunezcuVMyMjKi+w0ZMkSmTJkiept0+PDhMnDgQGndurVo0tavX79aEyyXsnknvhFUEgEEEEDAtwIu9bmhSdAef\/xxM3oWm6Dl5uZKy5Yto7cz9aGAzp07y9atWy9K0DQh07lqDRs2lH379kn\/\/v1NcB4\/ftyMvm3cuNEketVtLgWLb7+xnBgCCCCAgBMCLvW5oUnQCgsLze3IyNwyjdRhw4aZ0bM+ffpEA1f\/XUfLbrvtNvO306dPS\/v27aWsrMyU1yROHwyIbJqsTZs2zSR2JGhOfP+pJAIIIICATwVI0HzaMLWdVnl5ueg8tMWLF0tmZqZs2bJFxo0bJ5s2bTIPDOzfv988HDBv3jzRJzl1OQ3d5s+fb255aoK3cuVK8\/SnJmTXXXedeVhg4sSJ5hhNmjQhQQtgXHDKCCCAAALhESBBC2hbrl69WpYsWSIVFRXSuHFjmTFjhplztm3bNnniiSdkzZo1cubMGVmwYIH53x9\/\/LFJ2qZPny7Nmzc3c9T0AQJdS+3QoUNmRE1vk3bp0qVGEZeCJaBhwWkjgAACCIREwKU+NzS3OGNj7+TJk2aB2do2XeNMk7XYhwXqegzd36VgCcn3m2oggAACCARUwKU+N5QJWirjzqVgSaUrn4UAAggggEBVAZf6XBI0y\/h3KVgsqSiOAAIIIICAlYBLfS4JmlWocIvTko\/iCCCAAAIIeBYgQfNMxY4uBQutjQACCCCAQDoFXOpzGUGzjDSXgsWSiuIIIIAAAghYCbjU55KgWYUKtzgt+SiOAAIIIICAZwESNM9U7OhSsNDaCCCAAAIIpFPApT6XETTLSHMpWCypKI4AAggggICVgEt9LgmaVahwi9OSj+IIIIAAAgh4FiBB80zFji4FC62NAAIIIIBAOgVc6nMZQbOMNJeCxZKK4ggggAACCFgJuNTnkqBZhQq3OC35KI4AAggggIBnARI0z1Ts6FKw0NoIIIAAAgikU8ClPpcRNMtIcylYLKkojgACCCCAgJWAS30uCZpVqHCL05KP4ggggAACCHgWIEHzTMWOLgULrY0AAggggEA6BVzqcxlBs4w0l4LFkoriCCCAAAIIWAm41OeSoFmFCrc4LfkojgACCCCAgGcBEjTPVOzoUrDQ2ggggAACCKRTwKU+lxE0y0gLY7DMmjVLcnJyLGX8VZw6+as9ajob2ol2SpcAsZcu+bp9bhj73JoESNDqFhuX7B3GYKFOlkGRouK0U4qgLT+GdrIETFFx2ilF0JYfE8Z2IkHzGBQHDx6UFStWSGVlpQwYMEC6desmDRo0qLF0GIOFOnkMljTvRjuluQE8fjzt5BEqzbvRTmluAI8fH8Z2IkHz0PjHjh2TwYMHS8+ePaVv376ydOlSqV+\/vuTl5ZGgefDz8y5h\/FJTJz9H3P+dG+1EO6VLgNhLl3xiPpdbnDGOJSUlUlpaKvn5+eavp06dkl69esmyZcukRYsW1YrzBUhMICb7KLRTsoUTc3zaKTGOyT4K7ZRs4cQcn3ZKjGO6jkKCFiM\/adIk6dixo2RlZUX\/OnnyZHOrUxO16raHHnpIysrK0tV+fC4CCCCAAALOCNx+++1SXFzsRH1J0GKaOTs7W0aOHCndu3eP\/lWf7GnatKnov7EhgAACCCCAAAKpECBBi1F+\/PHHzehZbIKWm5srLVu2lKFDh6aiPfgMBBBAAAEEEEBASNBigqCwsFCaNWt2UTI2bNgwM3rWp08fwgUBBBBAAAEEEEiJAAlaDHN5ebnoPLTFixdLZmambNmyRcaNGyebNm2SjIyMlDQIH4IAAggggAACCJCgVYmB1atXy5IlS6SiokIaN24sM2bMkHbt2hEpCCCAAAIIIIBAygRI0GqgPnnypDRq1ChlDcEHIYAAAggggAACEQESNGIBAQQQQAABBBDwmQAJWpwNsmjRIjl8+LB5qGDEiBFxHsVfxd5++23585\/\/LDp62L59e3nwwQcDP4q4fv16+etf\/ypnz54VXT9H61SvXj1\/wcd5Nu+\/\/75om\/3oRz+K8wj+KHbgwAF55ZVX5MorrzQndO7cObn11lvl7rvv9scJWpyFvjpO30iiUyY6d+4sP\/jBDwL9ndJlh06fPh19\/Z22lf4nJycn0PXS+cd6nTh06JBZC1PfKFPbK\/4sQiJlRTds2GDqdObMGfnhD38oXbp0SdlnJ+qDtA7bt2+vNsY0DnVx+V27dknXrl2lf\/\/+gY7B6sxI0OKMpLfeekt2794ta9asMW8aCPq2bds2+elPfyq6rIiu+7Zw4UL58pe\/LM8++2xgq\/baa6\/Jb3\/7W1OnEydOyMyZM82FNwxLpmgSfc8998g3v\/lNWbBgQWDbSE9ck8xnnnnGrEGom154r732WunUqVOg6\/XZZ5\/JfffdZ+JNE845c+bIt771LROPQd10jq5u+go83bTzfOONN0T\/3rBhw0BWa9++ffLwww\/LL37xC3Pt06f5NZkeP358IOujJ718+XIpKiqSJ598Uj755BPRxFqv5bFLSAWhcu+9957oD7hRo0aZ64TOC49sY8eOFX09o66ysHPnTvNDaO3ataFK0kjQLKJ0x44d5j2dmswEfZsyZYrcdNNNMnz4cFOVDz\/8UAYNGmQuwEHdNEH7+te\/Hr0o6f\/XjkQ7yqBvEydOFE0ALly4IPPmzQt0dbSDf\/PNN+XXv\/51oOtR9eR1lF1\/xGnyqZuOuOuIgP5ICMOmo9I6Iv3zn\/9c7rjjjsBW6bnnnjM\/Rh955BFTh7\/\/\/e8mFv\/0pz8Ftk4DBw6UadOmRX\/krFy5Ul5++WX5\/e9\/H8g66Sur9LsUGdXUEekhQ4ZctMLC6NGjpV+\/fvLAAw8Eso7VnTQJmkVTaoL2\/PPPB34EQwn0VsyXvvQlc8s28gts8+bNpn5h2LQzmT59urmNNmHChEBXSS+y+gtZk2odIZw7d26g66O\/fDX+9Paf\/jDQW9FhWNZGO3wdFdRRzn\/84x9y8803S6tWrQLdVrEnr3Go7y7WEacgb\/oD+9133zXXB53+oPX629\/+Fn0ncxDrpgmNts0111xjTn\/\/\/v0mmQ7qD26tz549e6LTU\/T6p\/XTAZLIpj\/A9foR5JHPqrFGgmbx7QtTghZh0ERm\/vz58tJLL5lfXDqqFuTt888\/l8cee8wMgd9www2mTrHD5EGrm9ZDfynqBerf\/\/63uUAFPUHT2y+apH3ve98zPxB03qCuRxj0xaHvvPNOk5zpph2MjhTqCHUYbrHrlAGt3x\/+8AfzppUgb3pLXW9xfvDBB9K8eXMzMq2jZ1\/72tcCWy2Nsb59+5p40zmCulyUXtN1dLBJkyaBq1fVBG327Nlmbp3OfYxsOk1Hk229tRuWjQTNoiXDlqDpKIZ2\/vpLXzvIyGiaBVHai2rC+emnn8pHH31kFiDWL3rslzrtJ1iHE9D20WH9F154wazNp\/MzNEHTi5Xe6gzqpOatW7fK1VdfLR06dDAaeqHVkdug3o6JNKlOWu7Zs6dMnTrV\/EnnAumDD1q\/yAMRdWh+X+2qk7M3btwY+NEzRdUpD9omOpdJ4\/D11183IzE6t\/OKK67wlbvXk9Hz1zm3+t86T0vnq+p8ab3FHsQ6VU3Q9Ie2JtKx1\/J169aJjqIFfUQ3to1J0LxGfDX7hSlBO3XqlOhrrfSXVxientPm0gtTbJJZWVkp9957r5SVlQXyIqWjTNqZaCeim45ifPzxx+aXvv5C1ieZwrDpXC1NbPS2UxA7k0gbVDcnplu3bma0UEdzg7zpU4Fav6CPcp4\/f15uu+02Mzf1xhtvNE2iIzP6XdIFy9u0aRPkZhK9g6BTV44cOWJG0zRJC+JWNUHThLq4uFgKCgqi1dFro17zgz6FhQQtQRGqj2bn5+eHYg6a\/tr64osvor\/2lUiHxvVprSAuS6Hn3rt3b\/nZz35mkjIdYdKlHHRCvc6tC8OmCYyONAX9IYFf\/epX5gnAMWPGmIRML7pat6A\/naqd4YsvvmhuuehcII07nVCvc2eCnHhqJ6jzBKs+VRfU75Q+aasPbujotF7rdJ6W\/u\/YOVxBq9vTTz8tX\/nKV8ytW72G6xOqunRSUG+va4KmdwwiTw8fP37cPBCgdw90WRT98Z2VlWWuhfrEdFg2RtAsWlKzeL34huEpTv1FrE\/GREZn9FekDiFrpxLUuRj6hdbbSzrKpLeUdB0gXbMuLK\/u0g5SL1BBT9B0xEyTsb\/85S8m5vSC+5vf\/Eauv\/56i2+nP4rq\/Dqd+6Nzm\/S7pU\/Wffvb3\/bHycV5FjqPSdtHbyeFYdPJ57\/85S+j1wkdSdPOXp+EDOqmUzr0h45OH9BN74ro5Pkg3lrXJYV69OhhntiMnT+n83H1u6VPd+qPBp2WE6YnOLXdSNCC+g3kvD0L6BdcRyyCeHHyXMkQ7BhZ9DSo62nV1AQ6D1J\/8PDqOH8HqT4soDEYpnbSuNO5qZGRJ3+3QHxnF+bXMpKgxRcTlEIAAQQQQAABBJImQIKWNFoOjAACCCCAAALxCujc4cgcaB2JDvLczXgMSNDiUaMMAggggAACCCRVQN+A8Mc\/\/lG+8Y1vmIfYdA5nUJcTigeKBC0eNcoggAACCCCAgCcBnQunazZqcqXrHeoTmLGbPrW9atUq87T9\/fffb9ar1E3XO9P9I2skevqwEO1EghaixqQqCCCAAAII+E1AH77QF5nr4sa6KkDs65h27dplnq7X\/+gi6fruWl00WJ+k1WWE9MluLa9vrIi8L9Vv9UvW+ZCgJUuW4yKAAAIIIIBAVEDXBdS5ZPqO2simC8vqCFlkjTZd9kQTMV3iSRcV19E0XXD82WefFV0OKrKgsAusJGgutDJ1RAABBBBAIM0C1SVo+nYNfSuAjp5FNv2bvrpJR9e+853vSEZGhjz33HMyaNAgadWqVZprkbqPJ0FLnTWfhAACCCCAQOgE\/vnPf8qiRYvMGwsik\/jff\/99s5hx7KuXqiZoOpqmC4frorOahEU2fZPDlClT5D\/\/+Y95J68u8nzttdcG9j3K8TY4CVq8cpRDAAEEEEAAAXPbUt\/GoMnZU089ZV7Srrcx9R3BnTp1igpVTdD0oYDOnTubNx7EJmj6Fgd9qCDy2iYXl9hQNBI0vlwIIIAAAgggYCWgE\/mffPJJ0fdkvvPOOybBqvr0ZXW3OHXumY6W6UvrddM3Ouh7Q3X+mc49c3kjQXO59ak7AggggAACCRLQ5TIefPBB6d27txQWFl6ysKy+N\/j8+fPyk5\/8JPqJ+jedazZz5kzzt\/nz55tbnlre9Y0EzfUIoP4IIIAAAghYCmiSpctn6NIYr7\/+unz66acyffr0i96BPGfOHDNClpOTE\/00XSNtwYIFsmbNGvPCeh1103LNmze3PKPgFydBC34bUgMEEEAAAQTSJqAPCTz66KNSVFQkrVu3Nuehidrhw4clNzfX03npyJoma7Fz0TwVDPFOJGghblyqhgACCCCAQLIFdLL\/kSNH5Jprrrnoow4dOmQWpmWLT4AELT43SiGAAAIIIIAAAkkTIEFLGi0HRgABBBBAAAEE4hMgQYvPjVIIIIAAAggggEDSBEjQkkbLgRFAAAEEEEAAgfgESNDic6MUAggggAACCCCQNAEStKTRcmAEEEAAAQQQQCA+ARK0+NwohQACCCCAAAIIJE2ABC1ptBwYAQQQQAABBBCIT4AELT43SiGAAAIIIIAAAkkTIEFLGi0HRgABBBBAAAEE4hMgQYvPjVIIIIAAAggggEDSBEjQkkbLgRFAAAEEEEAAgfgESNDic6MUAggggAACCCCQNAEStKTRcmAEEEAAAQQQQCA+ARK0+NwohQACCCCAAAIIJE2ABC1ptBwYAQQQQAABBBCIT4AELT63aKnvf\/\/7UlZWZnmUuhW\/6qqr6lYgwXs3atQowUfkcAgggAACCFxe4Pbbb5fi4uLL7xiCPUjQLBuxRYsWkpWVZXmUuhXfsmVL3QrUce8dO3ZIx44dayy1devWOh4xsbv37t27zgfUJFq\/2InY4vn8RHyuHqNXr17RQ\/34xz+WRYsWJerQno+TzPrfcsstsnfvXs\/nEoQdqVMQWkmEdqKd\/CbgTIJ2+vRpKSkpkV27dknXrl2lf\/\/+Ut1I0IYNG2T79u1y7tw5ycnJqXaf2EbMzMw0X+wwbVUTtO9+97u+ql48Cerlks66VDDdCWpdzjVs+yYzOfRqVddzmDVrlrmWhGmjTsFozTC20\/z58+Wjjz4KRgNYnqUzCdrYsWPl2LFjkp2dLTt37pSlS5fK2rVrL0nA3nvvPTlw4ICMGjVK3n77bWncuHGtxPdkZsqtfkjQEphE\/e53v5OHH37YMrT8VZw6+as9ajob2qn2dornx0kiW54fJ4nUrNux6vrDQI8eljsHsVKadB45cqRueAHd24kEraKiQoYMGSKbNm2SjIwM01SjR4+Wfv36yQMPPFBt0+mo2O7du6VBgwa1Nm1uZqYM90OClubbjgGN\/\/hPu3v3+MsmqGT5jh3SqaC9ZgAAABkGSURBVJZb0Qn6mJQe5rJ1SuAPkVRVLNBJZw1x\/thjj8kLL7yQKsLo5yQzQfTSTluTPL3EVLSGGI8nOQ\/TnYNIgvrOO++QoKX8m5fED1y+fLmUlpZKXl5e9FNee+01+fDDD2X8+PE1Jmh79uyRevXq1XpmN910k1y4cKHGfYYPHx640SgvF6okNldSDp3wOvkgIS4vL5dOnTolxStdB71snXzgni6btHxuKhKStFTMxx9KjFfbOJGZzzuuukqeOnHCxw2YuFNzYgRt9uzZcubMmYvmgWzbtk0WLlwoRUVF1gnam2++mbgW4Ui+EkjmL3ZfVTRgJ0O7+LfB4hnp8W9twnFmYfq+6CoGJ0jQwhGYWouXX35ZPvvss4sStHXr1omOohUWFlolaIl6SCBMX6DwRI5IPPM+wlR\/v9aFdvFry1z8pLF\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\/v37zcPB8ybN0\/0SU5dTkO3+fPnm1uemuCtXLnSPP2pCdl1111nHhaYOHGiOUaTJk1I0AIYF5wyAggggEB4BEjQAtqWq1evliVLlkhFRYU0btxYZsyYYeacbdu2TZ544glZs2aNnDlzRhYsWGD+98cff2yStunTp0vz5s3NHDV9gEDXUjt06JAZUdPbpF26dKlRxKVgCWhYcNoIIIAAAiERcKnPDc0tztjYO3nypFlgtrZN1zjTZC32YYG6HkP3dylYQvL9phoIIIAAAgEVcKnPDWWClsq4cylYUunKZyGAAAIIIFBVwKU+lwTNMv5dChZLKoojgAACCCBgJeBSn0uCZhUq3OK05KM4AggggAACngVI0DxTsaNLwUJrI4AAAgggkE4Bl\/pcRtAsI82lYLGkojgCCCCAAAJWAi71uSRoVqHCLU5LPoojgAACCCDgWYAEzTMVO7oULLQ2AggggAAC6RRwqc9lBM0y0lwKFksqiiOAAAIIIGAl4FKfS4JmFSrc4rTkozgCCCCAAAKeBUjQPFOxo0vBQmsjgAACCCCQTgGX+lxG0CwjzaVgsaSiOAIIIIAAAlYCLvW5JGhWocItTks+iiOAAAIIIOBZgATNMxU7uhQstDYCCCCAAALpFHCpz2UEzTLSXAoWSyqKI4AAAgggYCXgUp9LgmYVKtzitOSjOAIIIIAAAp4FSNA8U7GjS8FCayOAAAIIIJBOAZf6XEbQLCPNpWCxpKI4AggggAACVgIu9bkkaFahwi1OSz6KI4AAAggg4FmABM0zFTu6FCy0NgIIIIAAAukUcKnPZQTNMtJcChZLKoojgAACCCBgJeBSn0uCZhUq3OK05KM4AggggAACngVI0DxTsaNLwUJrI4AAAgggkE4Bl\/pcRtAsI82lYLGkojgCCCCAAAJWAi71uSRoVqHCLU5LPoojgAACCCDgWYAEzTMVO7oULLQ2AggggAAC6RRwqc9lBM0y0lwKFksqiiOAAAIIIGAl4FKfS4JmFSrc4rTkozgCCCCAAAKeBUjQPFOxo0vBQmsjgAACCCCQTgGX+lxG0CwjzaVgsaSiOAIIIIAAAlYCLvW5JGhWocItTks+iiOAAAIIIOBZgATNMxU7uhQstDYCCCCAAALpFHCpz2UEzTLSXAoWSyqKI4AAAgggYCXgUp9LgmYVKuG8xTlr1izJycmxlPFXcerkr\/ao6WxoJ9opXQLEXrrk6\/a5JGh183J67zAGC3UKRkjTTrRTugSIvXTJ1+1zaae6efltb0bQqrTIwYMHZcWKFVJZWSkDBgyQbt26SYMGDWpsN74Afgvp6s+HdqKd0iVA7KVLvm6fSzvVzStde4exnWqyJEGLkTl27JgMHjxYevbsKX379pWlS5dK\/fr1JS8vjwQtXd\/GBH1uGL\/U1ClBwZHkw9BOSQZO0OFppwRBJvkwYWwnEjQPQVNSUiKlpaWSn59v9j516pT06tVLli1bJi1atKj2CGEMFurkIVh8sAvt5ING8HAKtJMHJB\/sQjv5oBE8nEIY24kEzUPDT5o0STp27ChZWVnRvSdPnmxudWqiVt320EMPSVlZmYejswsCCCCAAAII2AjcfvvtUlxcbHOIwJTlFmdMU2VnZ8vIkSOle\/fu0b\/qkz1NmzYV\/Tc2BBBAAAEEEEAgFQIkaDHKjz\/+uBk9i03QcnNzpWXLljJ06NBUtAefgQACCCCAAAIICAlaTBAUFhZKs2bNLkrGhg0bZkbP+vTpQ7gggAACCCCAAAIpESBBi2EuLy8XnYe2ePFiyczMlC1btsi4ceNk06ZNkpGRkZIG4UMQQAABBBBAAAEStCoxsHr1almyZIlUVFRI48aNZcaMGdKuXTsiBQEEEEAAAQQQSJkACVoN1CdPnpRGjRqlrCH4IAQQQAABBBBAICJAgkYsIIAAAggggAACPhMgQYuzQRYtWiSHDx82DxWMGDEizqP4q9jbb78tf\/7zn0VHD9u3by8PPvhg4EcR169fL3\/961\/l7NmzouvnaJ3q1avnL\/g4z+b9998XbbMf\/ehHcR7BH8UOHDggr7zyilx55ZXmhM6dOye33nqr3H333f44QYuz0FfH6RtJdMpE586d5Qc\/+EGgv1O67NDp06ejr7\/TttL\/5OTkBLpeOv9YrxOHDh0ya2HqG2Vqe8WfRUikrOiGDRtMnc6cOSM\/\/OEPpUuXLin77ER9kNZh+\/bt1caYxqEuLr9r1y7p2rWr9O\/fP9AxWJ0ZCVqckfTWW2\/J7t27Zc2aNeZNA0Hftm3bJj\/96U9FlxXRdd8WLlwoX\/7yl+XZZ58NbNVee+01+e1vf2vqdOLECZk5c6a58IZhyRRNou+55x755je\/KQsWLAhsG+mJa5L5zDPPmDUIddML77XXXiudOnUKdL0+++wzue+++0y8acI5Z84c+da3vmXiMaibztHVTV+Bp5t2nm+88Ybo3xs2bBjIau3bt08efvhh+cUvfmGuffo0vybT48ePD2R99KSXL18uRUVF8uSTT8onn3wimljrtTx2CakgVO69994T\/QE3atQoc53QeeGRbezYsaKvZ9RVFnbu3Gl+CK1duzZUSRoJmkWU7tixw7ynU5OZoG9TpkyRm266SYYPH26q8uGHH8qgQYPMBTiomyZoX\/\/616MXJf3\/2pFoRxn0beLEiaIJwIULF2TevHmBro528G+++ab8+te\/DnQ9qp68jrLrjzhNPnXTEXcdEdAfCWHYdFRaR6R\/\/vOfyx133BHYKj333HPmx+gjjzxi6vD3v\/\/dxOKf\/vSnwNZp4MCBMm3atOiPnJUrV8rLL78sv\/\/97wNZJ329k36XIqOaOiI9ZMiQi1ZYGD16tPTr108eeOCBQNaxupMmQbNoSk3Qnn\/++cCPYCiB3or50pe+ZG7ZRn6Bbd682dQvDJt2JtOnTze30SZMmBDoKulFVn8ha1KtI4Rz584NdH30l6\/Gn97+0x8Geis6DMvaaIevo4I6yvmPf\/xDbr75ZmnVqlWg2yr25DUO9d3FOuIU5E1\/YL\/77rvm+qDTH7Ref\/vb36LvZA5i3TSh0ba55pprzOnv37\/fJNNB\/cGt9dmzZ090eope\/7R+OkAS2fQHuF4\/gjzyWTXWSNAsvn1hStAiDJrIzJ8\/X1566SXzi0tH1YK8ff755\/LYY4+ZIfAbbrjB1Cl2mDxoddN66C9FvUD9+9\/\/NheooCdoevtFk7Tvfe975geCzhvU9QiDvjj0nXfeaZIz3bSD0ZFCHaEOwy12nTKg9fvDH\/5g3rQS5E1vqestzg8++ECaN29uRqZ19OxrX\/taYKulMda3b18TbzpHUJeL0mu6jg42adIkcPWqmqDNnj3bzK3TuY+RTafpaLKtt3bDspGgWbRk2BI0HcXQzl9\/6WsHGRlNsyBKe1FNOD\/99FP56KOPzALE+kWP\/VKn\/QTrcALaPjqs\/8ILL5i1+XR+hiZoerHSW51BndS8detWufrqq6VDhw5GQy+0OnIb1NsxkSbVScs9e\/aUqVOnmj\/pXCB98EHrF3kgog7N76tddXL2xo0bAz96pqg65UHbROcyaRy+\/vrrZiRG53ZeccUVvnL3ejJ6\/jrnVv9b52npfFWdL6232INYp6oJmv7Q1kQ69lq+bt060VG0oI\/oxrYxCZrXiK9mvzAlaKdOnRJ9rZX+8grD03PaXHphik0yKysr5d5775WysrJAXqR0lEk7E+1EdNNRjI8\/\/tj80tdfyPokUxg2nauliY3edgpiZxJpg+rmxHTr1s2MFupobpA3fSpQ6xf0Uc7z58\/LbbfdZuam3njjjaZJdGRGv0u6YHmbNm2C3EyidxB06sqRI0fMaJomaUHcqiZomlAXFxdLQUFBtDp6bdRrftCnsJCgJShC9dHs\/Pz8UMxB019bX3zxRfTXvhLp0Lg+rRXEZSn03Hv37i0\/+9nPTFKmI0y6lINOqNe5dWHYNIHRkaagPyTwq1\/9yjwBOGbMGJOQ6UVX6xb0p1O1M3zxxRfNLRedC6RxpxPqde5MkBNP7QR1nmDVp+qC+p3SJ231wQ0dndZrnc7T0v8dO4craHV7+umn5Stf+Yq5davXcH1CVZdOCurtdU3Q9I5B5Onh48ePmwcC9O6BLouiP76zsrLMtVCfmA7LxgiaRUtqFq8X3zA8xam\/iPXJmMjojP6K1CFk7VSCOhdDv9B6e0lHmfSWkq4DpGvWheXVXdpB6gUq6AmajphpMvaXv\/zFxJxecH\/zm9\/I9ddfb\/Ht9EdRnV+nc390bpN+t\/TJum9\/+9v+OLk4z0LnMWn76O2kMGw6+fyXv\/xl9DqhI2na2euTkEHddEqH\/tDR6QO66V0RnTwfxFvruqRQjx49zBObsfPndD6ufrf06U790aDTcsL0BKe2GwlaUL+BnLdnAf2C64hFEC9OnisZgh0ji54GdT2tmppA50HqDx5eHefvINWHBTQGw9ROGnc6NzUy8uTvFojv7ML8WkYStPhiglIIIIAAAggggEDSBEjQkkbLgRFAAAEEEEAgXgGdOxyZA60j0UGeuxmPAQlaPGqUQQABBBBAAIGkCugbEP74xz\/KN77xDfMQm87hDOpyQvFAkaDFo0YZBBBAAAEEEPAkoHPhdM1GTa50vUN9AjN206e2V61aZZ62v\/\/++816lbrpeme6f2SNRE8fFqKdSNBC1JhUBQEEEEAAAb8J6MMX+iJzXdxYVwWIfR3Trl27zNP1+h9dJF3fXauLBuuTtLqMkD7ZreX1jRWR96X6rX7JOh8StGTJclwEEEAAAQQQiArouoA6l0zfURvZdGFZHSGLrNGmy55oIqZLPOmi4jqapguOP\/vss6LLQUUWFHaBlQTNhVamjggggAACCKRZoLoETd+uoW8F0NGzyKZ\/01c36ejad77zHcnIyJDnnntOBg0aJK1atUpzLVL38SRoqbPmkxBAAAEEEAidwD\/\/+U9ZtGiReWNBZBL\/+++\/bxYzjn31UtUETUfTdOFwXXRWk7DIpm9ymDJlivznP\/8x7+TVRZ6vvfbawL5HOd4GJ0GLV45yCCCAAAIIIGBuW+rbGDQ5e+qpp8xL2vU2pr4juFOnTlGhqgmaPhTQuXNn88aD2ARN3+KgDxVEXtvk4hIbikaCxpcLAQQQQAABBKwEdCL\/k08+KfqezHfeecckWFWfvqzuFqfOPdPRMn1pvW76Rgd9b6jOP9O5Zy5vJGgutz51RwABBBBAIEECulzGgw8+KL1795bCwsJLFpbV9wafP39efvKTn0Q\/Uf+mc81mzpxp\/jZ\/\/nxzy1PLu76RoLkeAdQfAQQQQAABSwFNsnT5DF0a4\/XXX5dPP\/1Upk+fftE7kOfMmWNGyHJycqKfpmukLViwQNasWWNeWK+jblquefPmlmcU\/OIkaMFvQ2qAAAIIIIBA2gT0IYFHH31UioqKpHXr1uY8NFE7fPiw5ObmejovHVnTZC12LpqngiHeiQQtxI1L1RBAAAEEEEi2gE72P3LkiFxzzTUXfdShQ4fMwrRs8QmQoMXnRikEEEAAAQQQQCBpAiRoSaPlwAgggAACCCCAQHwCJGjxuVEKAQQQQAABBBBImgAJWtJoOTACCCCAAAIIIBCfAAlafG6UQgABBBBAAAEEkiZAgpY0Wg6MAAIIIIAAAgjEJ0CCFp8bpRBAAAEEEEAAgaQJkKAljZYDI4AAAggggAAC8QmQoMXnRikEEEAAAQQQQCBpAiRoSaPlwAgggAACCCCAQHwCJGjxuVEKAQQQQAABBBBImgAJWtJoOTACCCCAAAIIIBCfAAlafG6UQgABBBBAAAEEkiZAgpY0Wg6MAAIIIIAAAgjEJ0CCFp8bpRBAAAEEEEAAgaQJkKAljZYDI4AAAggggAAC8QmQoMXnRikEEEAAAQQQQCBpAiRoSaPlwAgggAACCCCAQHwCJGjxuVEKAQQQQAABBBBImsD\/B8UXGeAw1JZ4AAAAAElFTkSuQmCC","height":419,"width":559}}
%---
%[output:2148efe7]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":2,"type":"double","value":[["-0.199633804097405","0.114996362922158","0.099727310990734"],["-0.212743471372301","0.119489538233053","0.102662178221165"]]}}
%---
%[output:22350528]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":2,"type":"double","value":[["0.251044529478673"],["0.264720715464052"]]}}
%---
