clear all; close all; clc
%% Batch PreProcessing for the SNM lab
% Must be in outermost directory containing all study subject data
SubjectID = 'Func01';
SubjectDir = cd;
rootDir = fullfile(SubjectDir, SubjectID);

% Look for OTB files in subjects folder
filenames = dir(fullfile(rootDir,'*.otb'));
if isempty(filenames)
    filenames = dir(fullfile(rootDir,'*.otb+'));
end
[~,fIdx] = sort({filenames(:).date});
filenames = filenames(fIdx);

% Look for OTB files in subjects folder
c3dFiles = dir(fullfile(rootDir, '*.c3d'));
if isempty(c3dFiles)
    c3dFiles = '';
end
% [~,cIdx] = sort({c3dFiles(:).date});
% c3dFiles = c3dFiles(cIdx);
removeStatic = contains({c3dFiles.name},{'static','Static'});
staticC3D = c3dFiles(removeStatic,:);
c3dFiles(removeStatic,:) = [];

try
    staticTrial = fullfile(staticC3D.folder, staticC3D.name);
    static = Read_C3D(staticTrial);
    [staticR, staticMarker] = Get_StaticMarkers(static.Data.Markers);
    
    kkFolder = fullfile(rootDir, 'kk files');
    if exist(kkFolder,'dir') == 7
    elseif ~exist(kkFolder,'dir') && ~isempty(FileNamePathC3D)
        mkdir(kkFolder)
    else
    end
    
    saveStatic = fullfile(kkFolder,[SubjectID,'_static.mat']);
    save(saveStatic, 'staticR', 'staticMarker');
catch
end

%% key variables for preprocessing subfunctions
%%%% These must be properly filled out for the processing to work
%%%% automatically
Stickers = 3;   % number of arrays
AuxChans = 1;   % 13 for biodex setup, 1 for treadmill
MuscleList = {'TA', 'MG', 'Sol'}; % muscles tested, empty if multiple in skipped
Operator = 6; % who is preprocessing this data
TotalChans = Stickers * 64 + AuxChans;

% loop through files for preprocessing
for i = 1:length(filenames)
    FileNamePathOTB = fullfile(filenames(i).folder,filenames(i).name);
    try
        FileNamePathC3D = fullfile(c3dFiles(i).folder,c3dFiles(i).name);
    catch
        FileNamePathC3D = c3dFiles;
    end
    
    try
        disp(['Processing ', +newline, filenames(i).name, +newline])
        Format_SNMLab_Data(FileNamePathOTB, FileNamePathC3D, MuscleList, TotalChans, AuxChans, Operator, 2048)
    catch
        disp([filenames(i).name, ' Failed'])
    end
end
