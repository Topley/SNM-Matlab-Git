function [staticS] = staticC3dPro(rootDir, staticC3D)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
staticTrial = fullfile(staticC3D.folder, staticC3D.name);

kkFolder = fullfile(rootDir, 'kk files');
if exist(kkFolder,'dir') == 7
elseif ~exist(kkFolder,'dir') && ~isempty(FileNamePathC3D)
    mkdir(kkFolder)
else
end
[~,SubjectID] = fileparts(rootDir);

saveStatic = fullfile(kkFolder,[SubjectID,'_static.mat']);
if exist(saveStatic, 'file') ~= 2
    static = Read_C3D(staticTrial);
    [staticR, staticMarker] = Get_StaticMarkers(static.Data.Markers);
    staticS.staticR = staticR;
    staticS.staticMarker = staticMarker;
    save(saveStatic, 'staticR', 'staticMarker');
else
    try
        staticS = matfile(saveStatic);
    catch
        staticS = [];
    end
end

end

