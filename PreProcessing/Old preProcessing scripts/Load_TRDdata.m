
function [] = Load_TRDdata(FileNamePathOTB, FileNamePathC3D, MuscleList, ArrayIteration, TotalChans, matrixnames, Operator, fsamp)

%% DESCRIPTION

% UPDATES - This software has been updated to process files more
% efficiently and without the need to copy filenames into multiple lines.
% All file information will be extracted from the full file path that is
% the first input to the function "fullFilePath".
% This script has also been updated to pull in data from the c3d files used
% in the motion capture system to save the kinematic and kinetic data into
% the to_cluster file

%   Loads OT Biolab .otb file and unzips .sig file. Outputs channels with
%   possible saturation (y-value > 2000).  Plots total number of inputs;
%   multi-channel matrices are plotted in groups of 16. Prompts user for
%   bad channels to remove.  Saves data with bad channels removed.

%   Implemented file structure for data processing; threshold analysis
%   moved to 'plotMC'; reduced memory overhead by eliminating 'fopen'
%   command; remove saturated channels automatically

% UPDATED   Matt Topley
% AUTHOR    Laura Miller
%     Christopher Thompson

% DATE CREATED      6-Feb-2013
% DATE MODIFIED     22-May-2013
% DATE MODIFIED     06-Sept-2022

%% parse setup variables
muscle = MuscleList(ArrayIteration); % select muscle for filename
LastChan = ArrayIteration * 64; % calculate the last channel of the current array
FirstChan = LastChan - 63; % calculate the first channel of the current array
ArrayChannels = {FirstChan:LastChan}; % create cell array of the current array channels

%% Make trial directories if not created
[fileDir, otbTrialname] = fileparts(FileNamePathOTB);
[~, c3dTrialname] = fileparts(FileNamePathC3D);

toClusterFolder = fullfile(fileDir, 'to_cluster');
kkFolder = fullfile(fileDir, 'kk files');
trialFolder = fullfile(fileDir, otbTrialname);

if exist(trialFolder, 'dir') == 7
    OTBTrialPath = fullfile(trialFolder, [otbTrialname, '.otb']);
else
    mkdir (trialFolder); % create unzip directory
    copyfile(FileNamePathOTB, trialFolder);
    OTBTrialPath = fullfile(trialFolder, [otbTrialname, '.otb']);
    unzip(OTBTrialPath, trialFolder)
    delete(OTBTrialPath)
end

if exist(toClusterFolder,'dir') == 7
else
    mkdir(toClusterFolder)
end

if exist(kkFolder,'dir') == 7
    CreateKKFile = 0;
else
    mkdir(kkFolder)
    CreateKKFile = 1;
end

%% PreProcess EMG channels
[allbadchan, data] = Preprocess_OTB_EMG(OTBTrialPath, MuscleList, TotalChans, fsamp, matrixnames);

%% format EMG data for saving
if allbadchan == 'n'
    return
else
    allbadchanauto = 1:TotalChans;
    allbadchanauto(MuscleList{1}) = [];
    allbadchan = [allbadchan; allbadchanauto'];
    
    EMGall = data(MuscleList{1},:);
    
    goodchan=[];
    for i = 1:(TotalChans)
        if isempty(find(i==allbadchan))
            goodchan = [goodchan; i];
        end
    end
    
    % Quattro channels for Treadmill
    %%%%% LAST UPDATED 10/8/2022 %%%%%
    EMG = data(goodchan,:);
    otbEMG = data(TotalChans-1,:);
    trigger = data(TotalChans,:);
    
    %% save to_cluster with filename
    [~,subjectID] = fileparts(fileDir);
    saveid = strcat(subjectID,'_', otbTrialname(17:21),'_',muscle,'_',num2str(Operator), '.mat');
    try
        save(fullfile(toClusterFolder,saveid),'EMG','EMGall','fsamp','allbadchan', 'trigger', 'timeDelay', 'kkStop', '-v7.3','-nocompression');
    catch
        save(fullfile(toClusterFolder,saveid),'EMG','EMGall','fsamp','allbadchan', 'trigger', 'timeDelay', 'kkStop', '-v7.3','-nocompression');
    end
    
    %% read c3d file if kk file does not already exist
    if CreateKKFile == 1
        
        try
            c3dStruct =  Read_C3D(FileNamePathC3D);
            c3d_Data = c3dStruct.Data;
        catch
            disp('Cannot read trial c3d file')
        end
        
        FrameRate = c3dStruct.Headers.VideoFrameRate; % camera frame rate
        AVRatio = c3dStruct.Headers.AVRatio; % camera to analogs ratio
        TotalFrames = c3dStruct.Headers.LastFrame; % total camera frames
        fsampKinetic = AVRatio * FrameRate; % forceplate sampling rate
        c3dTime = [0:length(c3d_Data.EMG)-1]./960;
        
        try
            otbTime = [0:length(otbEMG)-1]./fsamp;
            try
                % synchronize forceplates and EMG data
                [fsKinetics, timeDiff] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, fsampKinetic, 0);
                % upsample forceplates to match EMG fsamp
                upsampPlates = resampler(cell2mat(struct2cell(c3d_Data.ForcePlates)'), fsKinetics, fsamp,0);
                
                % setup forceplate structure for easier data manipulation
                fpLabels = fieldnames(c3d_Data.ForcePlates);
                for ii = 1:size(fpLabels,1)
                    c3d_Data.ForcePlates.(fpLabels{ii}) = upsampPlates(:,ii);
                end
                
                 % upsample COP to match EMG fsamp
                upsampCOP = resampler(cell2mat(struct2cell(c3d_Data.COP)'), fsKinetics, fsamp,0);
                %setup COP structure for easier data manipulation
                copLabels = fieldnames(c3d_Data.COP);
                for jj = 1:size(copLabels,1)
                    c3d_Data.COP.(copLabels{jj}) = upsampCOP(:,jj);
                end
                
            catch
                disp('no kinetics in the trials c3d file')
            end
            
            try
                % synchronize forceplates and EMG data
                [fsVideo, timeDiffMarker] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, FrameRate, 0);
                % upsample marker trajectories to match EMG fsamp
                upsampMarkers = resampler(cell2mat(struct2cell(c3d_Data.Markers)'), fsVideo, fsamp,0);
                
                %setup Marker structure for easier data manipulation
                markerLabels = fieldnames(c3d_Data.Markers);
                for kk = 1:size(markerLabels,1)
                    c3d_Data.Markers.(markerLabels{kk}) = upsampMarkers(:,kk);
                end
                
            catch
                disp('no markers in the trials c3d file')
            end
            
            % creating timing variables to save in kk file
            delayLen = otbTime <= timeDiff;
            frontEnd = otbTime(delayLen);
            timeDelay = frontEnd(end);
            cutBack = length(upsampPlates);
            kkStop = cutBack./fsamp;
        catch
        end
        
        % save kk file
        subjID = otbTrialname(17:21);
        under = strfind(c3dTrialname, '_');
        saveC3D = [c3dTrialname(1:under(2)), subjID, '_kk.mat'];
        save(fullfile(kkFolder, saveC3D), 'c3d_Data', 'FrameRate', 'AVRatio', 'TotalFrames', 'trigger', 'timeDelay', 'fsKinetics', 'kkStop', '-v7.3','-nocompression');
    end
    
end
