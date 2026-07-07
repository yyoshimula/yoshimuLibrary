%[text] # Reading spacecraft shape data
%[text] 衛星形状や表面特性の読み込み
%[text] ## inputs
%[text] `satName` : satellite object file name, .obj file
%[text] ## outputs
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.normal`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.faces`: face indices, Nx3 matrix, 面を構成する座標のindex
%[text] `sat.area`: face area, m^2, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.uu`: x-axis of local frame, Nx3 matrix
%[text] `sat.uv`: y-axis of local frame, Nx3 matrix
%[text] `sat.qlb`: quaternion from body-fixed frame to local frame, Nx4 matrix
%[text] `sat.Ca`: coefficients for absorption, Nx1 vector, 吸収率 (MTL Red)
%[text] `sat.Cd`: coefficients for diffusion, Nx1 vector, 拡散反射率 (MTL Green)
%[text] `sat.Cs`: coefficients for specular reflection, Nx1 vector, 鏡面反射率 (MTL Blue)
%[text] `sat.F0`: reflectivity default, Nx1 vector
%[text] `sat.kappa`: Nx1 vector
%[text] `sat.nu`: default value for Ashikhmin-Shirley model, Nx1 vector
%[text] `sat.nv`: default value for Ashikhmin-Shirley model, Nx1 vector
%[text] `sat.mCT`: default value for Cook-Torrance model, Nx1 vector
%[text] `sat.fObs`: Nx1 vector
%[text] `sat.MOI`: moment of inertia, 3x3 matrix (default value)
%[text] `sat.m`: mass, kg (default value)
%[text] `sat.sunlitFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] `sat.force`: Nx3 matrix
%[text] `sat.torque`: Nx3 matrix
%[text] ## note
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20240729 major update, readObjで読み込むように変更 y.yoshimura
%[text] 2021020209  y.yoshimura
%[text] See also showSC.
function sat = readSC(satName)

sat = readObj(satName);

[sat.area, sat.pos] = calcAreaObj(sat);  % face areas and position
[sat.uu, sat.uv, sat.qlb] = calcLocalFrame(sat); % local frames and their quaternions

N = size(sat.area,1);
%[text] ### initialize
% reflectivity
sat.F0 = ones(N,1) .* 0.5; % default value
sat.kappa = zeros(N,1);
sat.Ca = zeros(N,1);
sat.Cd = ones(N,1) .* 0.5; % default value
sat.Cs = ones(N,1) .* 0.5; % default value

% Check for MTL file and parse if available
if isfield(sat, 'mtlFileName') && ~isempty(sat.mtlFileName)
    % The mtl file is expected to be in the same directory as the obj file
    [pathstr, ~, ~] = fileparts(satName);
    mtlFullPath = fullfile(pathstr, sat.mtlFileName);

    if exist(mtlFullPath, 'file')
        fid_mtl = fopen(mtlFullPath);

        % Initialize containers for material properties
        materials = containers.Map; % Name -> Struct('Ka',...,'Kd',...,'Ks',...)
        currentMatName = '';

        while 1
            tline = fgetl(fid_mtl);
            if ~ischar(tline), break, end

            tline = strtrim(tline);
            if isempty(tline) || tline(1) == '#', continue, end

            % Check for newmtl first
            if startsWith(tline, 'newmtl')
                currentMatName = sscanf(tline(8:end), '%s');
                materials(currentMatName) = struct('Ka', [0 0 0], 'Kd', [0.5 0.5 0.5], 'Ks', [0 0 0]);
                continue;
            end

            if ~isempty(currentMatName)
                matProps = materials(currentMatName);
                if startsWith(tline, 'Ka')
                    matProps.Ka = sscanf(tline(4:end), '%f')';
                elseif startsWith(tline, 'Kd')
                    matProps.Kd = sscanf(tline(4:end), '%f')';
                elseif startsWith(tline, 'Ks')
                    matProps.Ks = sscanf(tline(4:end), '%f')';
                end
                materials(currentMatName) = matProps;
            end
        end
        fclose(fid_mtl);

        % Assign values to sat.Ca, sat.Cd, sat.Cs
        % Map RGB of Kd to Ca, Cd, Cs respectively (Nx1 vectors)
        sat.Ca = zeros(N, 1);
        sat.Cd = zeros(N, 1);
        sat.Cs = zeros(N, 1);

        for k = 1:length(sat.materialNames)
            matName = sat.materialNames{k};
            if isKey(materials, matName)
                props = materials(matName);
                % Find faces using this material
                faceIndices = find(sat.materialIndices == k);

                % Assign properties: R -> Ca, G -> Cd, B -> Cs
                % props.Kd is [R, G, B]
                sat.Ca(faceIndices, 1) = props.Kd(1);
                sat.Cd(faceIndices, 1) = props.Kd(2);
                sat.Cs(faceIndices, 1) = props.Kd(3);
            end
        end
        disp(['MTL file parsed and applied RGB->Ca,Cd,Cs: ' mtlFullPath]);
    else
        disp(['Warning: MTL file specified but not found: ' mtlFullPath]);
    end
end

sat.nu = ones(N,1) .* 800; % default value for Ashikhmin-Shirley model
sat.nv = ones(N,1) .* 800; % default value for Ashikhmin-Shirley model
sat.mCT = ones(N,1) .* 0.05; % default value for Cook-Torrance model
sat.fObs = zeros(N,1);

% physical (default values)
sat.MOI = diag([10, 15, 20]);
sat.m = 100;

sat.sunlitFlag = ones(N,1);

sat.force = zeros(N,3);
sat.torque = zeros(N,3);

disp('number of vertices:')
disp(size(sat.vertices, 1)) % nVertices
disp('number of faces:')
disp(size(sat.faces, 1)) % nFaces

end

%[appendix]{"version":"1.0"}
%---
