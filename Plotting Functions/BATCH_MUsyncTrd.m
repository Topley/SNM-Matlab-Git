rootDir = cd;
% subC3ds = dir('MartinTest\*.c3d');
% [~,idx] = sort({subC3ds.date});
% subC3ds = subC3ds(idx);

matFiles = dir('MartinTest\decomposed\*5.mat');
clusterFileVariables = {'fsamp','CrossCorrEMG', 'AnalogData', 'AnalogNames'};

for i = 1%:size(matFiles, 1)
    sVars = load(fullfile(matFiles(i).folder,matFiles(i).name), clusterFileVariables{:});
    decompFile = [matFiles(i).name(1:end-4), '_v23decomposed.mat'];
    filePulses = load(fullfile(matFiles(i).folder, decompFile), 'MUPulses');
    emgChanC3d = contains(sVars.AnalogNames(:), 'EMG');
    c3dEMG = sVars.AnalogData(emgChanC3d,:);
    otbEMG = sVars.CrossCorrEMG;
    try
        [trueFs,diffT] = emg_sync(otbEMG, double(c3dEMG), sVars.fsamp, 960, 1)
        
    catch
        disp('C3D file not found')
    end
    
end
fsamp = 2048;

rFxChanC3d = contains(sVars.AnalogNames(:), '.rFx');
rFx = double(sVars.AnalogData(rFxChanC3d,:));
rFx = rFx-round(mean(1:500));

rFzChanC3d = contains(sVars.AnalogNames(:), '.rFz');
rFz = double(sVars.AnalogData(rFzChanC3d,:));
rFz = rFz-round(mean(1:500));
timeRFz = [1:length(rFz)]./trueFs;

% timeIdx = find(timeRFz >= diffT);
% adjustedTime = timeRFz(adjustedTime);
MUFiring = filePulses.MUPulses;
for i = 1:length(MUFiring)
MUFiring{i}(find(MUFiring{i} <= diffT*trueFs)) = [];

MUFiring{i} = MUFiring{i}-MUFiring{i}(1);
end 


figure
plot(timeRFz, rFz)
hold on
for i =1:length(MUFiring)
    plot(MUFiring{i}(2:end)./fsamp, (1./diff(MUFiring{i}))*fsamp, 'o')
end 