%[text] # Reading .obj file 
%[text] wavefront obj fileの読み込み
%[text] `sat`: 構造体satに各種データを格納
%[text] `sat.vertices`: vertex position, (x, y, z), m, Nx3 matrix, 面を構成する点の座標
%[text] `sat.faces`: face indices, Nx3 (or Nx4) matrix, 面を構成する座標のindex
%[text] `sat.area`: face area, ${\\rm m^2}$, Nx1 vector, 面の面積
%[text] `sat.pos`: center of face, m, Nx3 matrix, 面の中心（平均）座標
%[text] `sat.normal`: normal vector, Nx3 matrix, 法線方向単位ベクトル, 外向きが正
%[text] `sat.shadowFlag`: self shadow flag, (default)1: not shadowoed, 0: shadowed
%[text] ## note
%[text] NA
%[text] ## references
%[text] modified following
%[text] obj = readObj(fname)
%[text] This function parses wavefront object data It reads the mesh vertices, texture coordinates, normal coordinates and face definitions(grouped by number of vertices) in a .obj file
%[text] INPUT: fname - wavefront object file full path
%[text] OUTPUT: obj.v - mesh vertices : obj.vt - texture coordinates : obj.vn - normal coordinates : obj.f - face definition assuming faces are made of of 3 vertices
%[text] Bernard Abayowa, Tec^Edge 11/8/07 
%[text] ## revisions
%[text] 20240729  modified by Yasuhiro Yoshimura
%[text] See also showSC, readSC.
function obj = readObj(fname)

%[text] ## pre-allocation
v = [];
vt = [];
vn = [];
f.v = [];
f.vt = [];
f.vn = [];

fid = fopen(fname);

%[text] ## parse .obj file
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end  % exit at end of file

    ln = sscanf(tline,'%s',1); % line type
    %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            v = [v
                sscanf(tline(2:end),'%f')'];

        case 'vt'  % texture coordinate
            vt = [vt
                sscanf(tline(3:end),'%f')'];

        case 'vn'  % normal coordinate
            vn = [vn
                sscanf(tline(3:end),'%f')'];

        case 'f'   % face definition
            fv = []; fvt = []; fvn = [];
            str = textscan(tline(2:end),'%s');
            str = str{1};

            nf = length(strfind(str{1},'/')); % number of fields with this face vertices

            [tok, str] = strtok(str,'//');     % vertex only
            for k = 1:length(tok)
                fv = [fv str2num(tok{k})];
            end

            if (nf > 0)
                [tok, str] = strtok(str,'//');   % add texture coordinates
                for k = 1:length(tok)
                    fvt = [fvt str2num(tok{k})];
                end
            end

            if (nf > 1)
                [tok, str] = strtok(str,'//');   % add normal coordinates
                for k = 1:length(tok)
                    fvn = [fvn str2num(tok{k})];
                end
            end

            if length(fv)<4
                length_NaN_vector = 4-length(fv);
                fv = [fv, NaN(1,length_NaN_vector)];
                fvt = [fvt, NaN(1,length_NaN_vector)];
                fvn = [fvn, NaN(1,length_NaN_vector)];
            end
            
            f.v = [f.v; fv]; f.vt = [f.vt; fvt]; f.vn = [f.vn; fvn];
    end
end
fclose(fid);

%[text] ## output
normalVec = vn(f.vn(:,1),:);

obj.vertices = v; % coordinates of vertices
obj.normal = normalVec; % normal vector of faces
obj.faces = f.v; % face indices

% obj.vt = vt; % for texture

end


%[appendix]{"version":"1.0"}
%---
