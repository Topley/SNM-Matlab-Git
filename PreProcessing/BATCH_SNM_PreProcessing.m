clear all
close all
clc
%% Batch PreProcessing for the SNM lab
    % Must be in outermost directory containing all study subject data
SubjectID = 'Func01';
SubjectDir = cd;
rootDir = fullfile(SubjectDir, SubjectID);

[filenames, c3dFiles, staticC3D] = SortTrials(rootDir);
staticMarkers = staticC3dPro(rootDir, staticC3D);

%% KEY VARIABLES FOR PROCESSING 
Stickers = 3;   % number of arrays
AuxChans = 1;   % 13 for biodex setup, 1 for treadmill
MuscleList = {'TA', 'MG', 'Sol'}; % muscles tested, empty if multiple in skipped
TotalChans = Stickers * 64 + AuxChans;

%% loop through files for preprocessing

for i = 1:length(filenames)
    FileNamePathOTB = fullfile(filenames(i).folder,filenames(i).name);
    FileNamePathC3D = fullfile(c3dFiles(i).folder,c3dFiles(i).name);
    try
        disp(['Processing ', +newline, filenames(i).name, +newline])
        Format_SNMLab_Data(MuscleList, TotalChans, AuxChans, FileNamePathOTB, FileNamePathC3D, staticMarkers)
    catch
        disp([filenames(i).name, ' Failed'])
    end
end

%% Treadmill and Coco file identifiers 
% 0 = max Forward  ||  max PF
% 1 = max Backward  ||  max DF
% 2 = Forward Lean Holds  ||  12 = 10% PF Hold  ||  22 = 20% PF Hold
% 3 = Backward Lean Holds  ||  13 = 10% DF Hold  ||  23 = 20% DF Hold
% 4 = Forward Lean Ramps  ||  14 = 10% PF Ramps  ||  24 = 20% PF Ramps
% 5 = Backward Lean Ramps  ||  15 = 10% DF Ramps  ||  25 = 20% DF Ramps
% 6 = PF EMG Holds  ||  16 = 10% Coco Holds  ||  26 = 20% Coco Holds
% 7 = DF EMG Holds  ||  17 = 10% Coco Ramps  ||  27 = 20% Coco Ramps
% 8 = Forward Ankle Angle ||  18 = 10% Coco Inverted Ramps  ||  28 = 20% Coco Inverted Ramps
% 9 = Backward Ankle Angle
% 10 = Quiet Stance  
% 11 = LoS
