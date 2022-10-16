
function [] = Format_TRDdata(FileNamePathOTB, FileNamePathC3D, MuscleList, ArrayIteration, TotalChans, AuxChans, Operator, fsamp)

%% DESCRIPTION

% UPDATES 10/9/2022
% This function has been modified to process EMG data from the OTB
% amplifier and Kinetic and Kinematic data from Vicon/Bertec.
% Folder and file creation has been modified for faster processing
% This function will extract EMG data from the otb amplifier in the old otb
% format and the new otb+ format. C3d is the standard file format for 3d
% data so the c3d_reader should be universal for major motion capture
% systems

%   Loads OT Biolab .otb file and unzips .sig file. This is also updated for the new file format
%   Outputs channels with possible saturation (y-value > 2000).  Plots total number of inputs;
%   multi-channel matrices are plotted in groups of 16. Prompts user for
%   bad channels to remove.  Saves data with bad channels removed.

%   Implemented file structure for data processing; threshold analysis
%   moved to 'Preprocess_OTB_EMG'; reduced memory overhead by eliminating 'fopen'
%   command; remove saturated channels automatically

%   If c3d data is present, it will be extracted, synched with EMG, and
%   resampled to 2048Hz. All kinetic and kinematic data will be saved in a
%   'kk file' using the same naming structure as the to_cluster without the
%   muscle tag

% UPDATED   Matt Topley
% AUTHOR    Laura Miller
%     Christopher Thompson

% DATE CREATED      6-Feb-2013
% DATE MODIFIED     22-May-2013
% DATE MODIFIED     06-Sept-2022

%% parse setup variables
Muscle = MuscleList{ArrayIteration}; % select muscle for filename
if contains(Muscle, 'LG')
    return
end

LastChan = ArrayIteration * 64; % calculate the last channel of the current array
FirstChan = LastChan - 63; % calculate the first channel of the current array
ArrayChannels = {FirstChan:LastChan}; % create cell array of the current array channels

%% Make trial directories if not created

% get OTB file information
[fileDir, otbTrialname] = fileparts(FileNamePathOTB);
[~, c3dTrialname] = fileparts(FileNamePathC3D);

% create directory names
toClusterFolder = fullfile(fileDir, 'to_cluster');
kkFolder = fullfile(fileDir, 'kk files');
trialFolder = fullfile(fileDir, otbTrialname);

% check if OTB file has a specific file folder. Create one if it does not
if exist(trialFolder, 'dir') == 7
else
    mkdir (trialFolder);
    copyfile(FileNamePathOTB, trialFolder);
    OTBTrialPath = fullfile(trialFolder, [otbTrialname, '.otb']);
    
    % Check which OTB version file is saved as. The new OTB files need to be opened differently
    if exist(OTBTrialPath)
        unzip(OTBTrialPath, trialFolder)
        delete(OTBTrialPath)
    else
        OTBTrialPath = fullfile(trialFolder, [otbTrialname, '.otb+']);
        untar(OTBTrialPath, trialFolder)
        delete(OTBTrialPath)
    end
    
end

% make to_cluster folder if not existent
if exist(toClusterFolder,'dir') == 7
else
    mkdir(toClusterFolder)
end

% make kk folder if not existent and c3d files exist
if exist(kkFolder,'dir') == 7
elseif ~exist(kkFolder,'dir') && ~contains(FileNamePathC3D, 'None')
    mkdir(kkFolder)
else
end

%% create to_cluster and kk filenames for saving
[~,subjectID] = fileparts(fileDir);
saveid = strcat(subjectID, '_', otbTrialname(end-4:end), '_', Muscle, '_', num2str(Operator), '.mat');

% save kk file
save_KKFile = erase([saveid(1:end-4), '_kk.mat'], ['_', Muscle]);
if exist(fullfile(kkFolder, save_KKFile))
    CreateKKFile = 0;
else
    CreateKKFile = 1;
end

%% PreProcess EMG channels
[allbadchan, data] = Preprocess_OTB_EMG(trialFolder, Muscle, ArrayChannels, TotalChans, fsamp);

%% format EMG data for saving
if allbadchan == 'n'
    return
else
    allbadchanauto = 1:TotalChans;
    allbadchanauto(ArrayChannels{:}) = [];
    allbadchan = [allbadchan; allbadchanauto'];
    
    EMGall = data(ArrayChannels{:},:);
    
    goodchan=[];
    for i = 1:(TotalChans)
        if isempty(find(i==allbadchan))
            goodchan = [goodchan; i];
        end
    end
    
    %%%%% LAST UPDATED 10/8/2022 %%%%%
    % this is always saved at this index
    EMG = data(goodchan,:);
    
    if AuxChans == 13
        
        %%%% Quattro Biodex Standrard Setup %%%%
        % Make sure BNC cables are properly connected when collecting for
        % these to be accurate
        
        Torque = data(totalchan-12,:);% 257 correct 3
        Velocity = data(totalchan-11,:);% 258 incorrect 4 60? possibly 180-120(JTAngle)
        Position = data(totalchan-10,:);% 259 incorrect 5 23960 ??
        Trigger = data(totalchan-9,:);% 260 correct 6
        EMGTorque = data(totalchan-8,:);% 261 EMG 7
        JR3Fx = data(totalchan-7,:);% 262 correct 8
        JR3Fy = data(totalchan-6,:);% 263 correct 9
        JR3Fz = data(totalchan-5,:);% 264 correct 10
        JR3Mx = data(totalchan-4,:);% 265 correct 11
        JR3My = data(totalchan-3,:);% 266 correct 12
        JR3Mz = data(totalchan-2,:);% 267 correct 13
        TargetEMG = data(totalchan-1,:);% 268 no feedback 14
        TargetTorque = data(totalchan,:); % 269 no feedback 15
        
        try
            save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'allbadchan',  'Torque', 'Velocity', 'Position', 'Trigger', 'EMGTorque', ...
                'JR3Fx','JR3Fy','JR3Fz','JR3Mx','JR3My','JR3Mz', 'TargetEMG', 'TargetTorque', '-v7.3','-nocompression');
            saved = [saveID, 'Complete']
            
        catch
            disp(['Unable to save Biodex setup for ', saveid])
        end
    end
    
    % If not using Biodex, the file will be saved as this
    otbEMG = data(TotalChans-1,:); % emg channel plugged into from analog out to Aux1 for synchronizing
    try
        save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'otbEMG', 'allbadchan', '-v7.3', '-nocompression');
        saved = [saveID, 'Complete']
    catch
        disp(['Unable to save ', saveid])
    end
    
    %% read c3d file 
    % if kk file does not already exist and there is c3d data available
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
        fsampPlates = AVRatio * FrameRate; % forceplate sampling rate
        c3dTime = [0:length(c3d_Data.EMG)-1]./fsampPlates;
        
        otbTime = [0:length(otbEMG)-1]./fsamp;
        
        try
            try
                
                % synchronize forceplates and EMG data
                [fsampK, timeDiff] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, fsampPlates, 0);
                
                % upsample forceplates to match EMG fsamp
                upsampPlates = resampler(cell2mat(struct2cell(c3d_Data.ForcePlates)'), fsampK, fsamp, 0);
                
                % setup forceplate structure for easier data manipulation
                fpLabels = fieldnames(c3d_Data.ForcePlates);
                for ii = 1:size(fpLabels,1)
                    c3d_Data.ForcePlates.(fpLabels{ii}) = upsampPlates(:,ii);
                end
                ForcePlates = c3d_Data.ForcePlates;
                
                % upsample COP to match EMG fsamp
                upsampCOP = resampler(cell2mat(struct2cell(c3d_Data.COP)'), fsampK, fsamp, 0);
                
                %setup COP structure for easier data manipulation
                copLabels = fieldnames(c3d_Data.COP);
                for jj = 1:size(copLabels,1)
                    c3d_Data.COP.(copLabels{jj}) = upsampCOP(:,jj);
                end
                COP = c3d_Data.COP;
                
            catch
                disp('no kinetics in the trials c3d file')
            end
            
            try
                % synchronize forceplates and EMG data
                
                [fsampMars, timeDiffMarker] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, FrameRate, 0);
                % upsample marker trajectories to match EMG fsamp
                upsampMarkers = resampler(cell2mat(struct2cell(c3d_Data.Markers)'), fsampMars, fsamp, 0);
                
                %setup Marker structure for easier data manipulation
                markerLabels = fieldnames(c3d_Data.Markers);
                for kk = 1:size(markerLabels,1)
                    c3d_Data.Markers.(markerLabels{kk}) = upsampMarkers(:,kk);
                end
                Markers = c3d_Data.Markers;
                
            catch
                disp('no markers in the trials c3d file')
            end
            
            % creating timing variables to save in kk file
            delayLen = otbTime <= timeDiff;
            frontEnd = otbTime(delayLen);
            TimeDelay = frontEnd(end);
            cutBack = length(upsampPlates);
            kkStop = cutBack./fsamp;
        catch
        end
        
        try
            save(fullfile(kkFolder, save_KKFile), 'c3d_Data', 'Markers', 'fsampMars', 'COP', 'ForcePlates', 'fsampK', 'FrameRate', 'AVRatio', 'TotalFrames', 'TimeDelay', 'kkStop', '-v7.3','-nocompression');
        catch
            save(fullfile(kkFolder, save_KKFile), 'c3d_Data', 'COP', 'ForcePlates','fsampK', 'AVRatio', 'TotalFrames', 'TimeDelay', 'kkStop', '-v7.3','-nocompression');
        end
    else
    end
    
end
