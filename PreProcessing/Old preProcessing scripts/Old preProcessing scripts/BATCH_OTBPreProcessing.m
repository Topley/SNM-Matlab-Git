clear all
close all
clc
%% Batch PreProcessing for the SNM lab

% Choose the outermost directory with all study subject data
EnterSubject = 'TRDTest04';
SubjectDir = cd;
rootDir = fullfile(SubjectDir, EnterSubject);

% Look for OTB files in subjects folder
filenames = dir(fullfile(rootDir,'*.otb'));
if isempty(filenames)
    % If the new software was used, find the correct fiels
    filenames = dir(fullfile(rootDir,'*.otb+'));
end

% Look for OTB files in subjects folder
c3dFiles = dir(fullfile(rootDir, '*.c3d'));
if isempty(c3dFiles)
    % If no c3d files exist in subjects folder then ignore kk file creation
    c3dFiles = 'None';
end

% removes the static trial used to build realtime feedback model of
% kinematics 
removeStatic = contains({c3dFiles.name},{'static','Static'});
c3dFiles(removeStatic,:) = [];

%% key variables for preprocessing subfunctions
%%%% These must be properly filled out for the processing to work
%%%% automatically

Stickers = 4;   % number of arrays
AuxChans = 1;   % 13 for biodex setup, 0 for treadmill
MuscleList = {'TA', 'MG', 'LG', 'Sol'}; % muscles tested
Operator = 5; % who is preprocessing this data
TotalChans = Stickers * 64 + AuxChans;

% loop through files for preprocessing
    for i = 1:length(filenames)
            FileNamePathOTB = fullfile(filenames(i).folder,filenames(i).name);
            
            % making sure preprocessing runs if c3ds are missing
            try
            FileNamePathC3D = fullfile(c3dFiles(i).folder,c3dFiles(i).name);
            catch
                FileNamePathC3D = c3dFiles; 
            end 
            
            try
                for j = 1:size(MuscleList,2)
                    disp([+newline, 'Processing ', +newline, FileNamePathOTB,'    &     ', FileNamePathC3D])
                ArrayIteration = j;
                Format_TRDdata(FileNamePathOTB, FileNamePathC3D, MuscleList, ArrayIteration, TotalChans, AuxChans, Operator, 2048)
                end 
            catch
                disp([filenames(i).name, ' Failed'])
            end
    end
