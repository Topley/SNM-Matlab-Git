% addpath('F:\Toolkit\Mirror\TRD Hub\Tribal Realtime\627 Biomechanical methods');
%  addpath('C:\Users\tuc43377\builds\c3d executable\ReadC3D_Standalone2');
clearvars
tic
%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
subC3ds = dir('MartinTest\*.c3d');
[~,idx] = sort({subC3ds.date});
subC3ds = subC3ds(idx);

matFiles = dir('MartinTest\decomposed\*5.mat');

for i = 1:size(subC3ds, 1)
    
    idxT = i*3;
    loopMatFile = matFiles(idxT-2:idxT,:);
    loopC3dFile = subC3ds(i,:);
    
    try
        [AnalogData, AnalogNames] = loopC3d(loopC3dFile);
        for j = 1:size(loopMatFile,1)
            try
            save(fullfile(loopMatFile(j).folder, loopMatFile(j).name), 'AnalogData', 'AnalogNames', '-append')
            catch
                disp(['Failed to process ', loopMatFile(j).name])
            end
        end 
    catch
        disp('C3D file not found')
    end
    
end
toc

function [AnalogData, AnalogNames] = loopC3d(fullFilename)

%[folderLoc,~] = fileparts(fullFilename);
e=actxserver('Readc3dLV.Application');
vipath= 'C:\SNM Matlab SDK\Readc3dLV.exe\ReadC3D_Standalone.vi';
vi=invoke(e,'GetVIReference',vipath);
% test = e.GetVIReference(vipath, [true])%, [0x08])

% cd(folderLoc);
vi.SetControlValue('selected path', fullfile(fullFilename.folder, fullFilename.name));
vi.Run
AnalogData = vi.GetControlValue('Analog Data');
AnalogNames = vi.GetControlValue('Analog Info');
delete(vi)

AnalogNames(13:14,:) = [];
AnalogData(13:14,:) = [];

end

