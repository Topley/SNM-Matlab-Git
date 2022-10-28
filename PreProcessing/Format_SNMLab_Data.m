
function [] = Format_SNMLab_data(FileNamePathOTB, FileNamePathC3D, MuscleList, TotalChans, AuxChans, Operator, fsamp)

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

ArrayNumber = size(MuscleList,2);

%% Make trial directories if not created
% get OTB file information
[fileDir, otbTrialname] = fileparts(FileNamePathOTB);

% get c3d information or insert ''
if ~isempty(FileNamePathC3D)
    [~, c3dTrialname] = fileparts(FileNamePathC3D);
end

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
elseif ~exist(kkFolder,'dir') && ~isempty(FileNamePathC3D)
    mkdir(kkFolder)
else
end

%% Extract data from OTB sig file
%%%% opening the sig file is different with OTB+ files. The script is updated to open old and new OTB files %%%%

sigFile = dir(fullfile(trialFolder,'*.sig'));
sigFile = fullfile(sigFile.folder, sigFile.name);
f = fopen(sigFile);

%%% extract data from Quattro
ChannelMatrix = fread(f,[TotalChans+8,inf],'short');
fclose all;
%% Loop through and PreProcess EMG channels
 for ij = 1:ArrayNumber
    
     Muscle = MuscleList{4} % select muscle for filename
     
     % skip preprocessing muscle if array is empty
    if isempty(Muscle)
        continue
    else
        LastChan = ij * 64; % calculate the last channel of the current array
        FirstChan = LastChan - 63; % calculate the first channel of the current array
        ArrayChannels = {FirstChan:LastChan}; % create cell array of the current array channels
        
        % create to_cluster and kk filenames for saving
        [~,subjectID] = fileparts(fileDir);
        if ~contains(FileNamePathOTB, '_') && contains(FileNamePathOTB, '.otb+')
            newName = [otbTrialname,'00'];
            saveid = strcat(subjectID, '_', newName(end-4:end), '_', Muscle, '_', num2str(Operator), '.mat');
        elseif contains(otbTrialname, '_')
            newName = erase(otbTrialname,'_');
            saveid = strcat(subjectID, '_', newName(end-4:end), '_', Muscle, '_', num2str(Operator), '.mat');
        else
            saveid = strcat(subjectID, '_', otbTrialname(end-4:end), '_', Muscle, '_', num2str(Operator), '.mat');
        end
        
        save_KKFile = replace(saveid, Muscle, 'kk'); % save kk file
        
        % preprocess OTB data
        [allbadchan, data] = Preprocess_OTB_EMG(ChannelMatrix, Muscle, ArrayChannels, TotalChans, fsamp);
        
        %% format EMG data for saving
        if allbadchan == 'n'
            return
        else
            allbadchanauto = 1:TotalChans;
            allbadchanauto(ArrayChannels{:}) = [];
            allbadchan = [allbadchan; allbadchanauto'];
            
            goodchan=[];
            for ik = 1:(TotalChans)
                if isempty(find(ik == allbadchan))
                    goodchan = [goodchan; ik];
                end
            end
            
            %%%%% LAST UPDATED 10/8/2022 %%%%%
            % this is always saved at this index
            EMGall = data(ArrayChannels{:},:);
            EMG = data(goodchan,:);
            
            if AuxChans == 13
                
                %%%% Quattro Biodex Standrard Setup %%%%
                % Make sure BNC cables are properly connected when collecting for
                % these to be accurate
                
                Torque = data(TotalChans-12,:);% 257 correct 3
                Velocity = data(TotalChans-11,:);% 258 incorrect 4 60? possibly 180-120(JTAngle)
                Position = data(TotalChans-10,:);% 259 incorrect 5 23960 ??
                Trigger = data(TotalChans-9,:);% 260 correct 6
                EMGTorque = data(TotalChans-8,:);% 261 EMG 7
                JR3Fx = data(TotalChans-7,:);% 262 correct 8
                JR3Fy = data(TotalChans-6,:);% 263 correct 9
                JR3Fz = data(TotalChans-5,:);% 264 correct 10
                JR3Mx = data(TotalChans-4,:);% 265 correct 11
                JR3My = data(TotalChans-3,:);% 266 correct 12
                JR3Mz = data(TotalChans-2,:);% 267 correct 13
                TargetEMG = data(TotalChans-1,:);% 268 no feedback 14
                TargetTorque = data(TotalChans,:); % 269 no feedback 15
                
                try
                    save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'allbadchan',  'Torque', 'Velocity', 'Position', 'Trigger', 'EMGTorque', ...
                        'JR3Fx','JR3Fy','JR3Fz','JR3Mx','JR3My','JR3Mz', 'TargetEMG', 'TargetTorque', '-v7.3','-nocompression');
                catch
                    disp(['Unable to save Biodex setup for ', saveid])
                end
            end
            
            % If not using Biodex, the file will be saved as this
            otbEMG = data(TotalChans,:); 
            otbEMG = ChannelMatrix(TotalChans,:);% emg channel plugged into from analog out to Aux1 for synchronizing
            trialLength = length(otbEMG);
            try
                save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'otbEMG', 'allbadchan', '-v7.3', '-nocompression');
            catch
                save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'allbadchan', '-v7.3', '-nocompression');
            end
           
            %% read c3d file
            % if kk file does not already exist and there is c3d data available
            if ij == 4
                try
                    c3dStruct =  Read_C3D(FileNamePathC3D);
                    c3d_Data = c3dStruct.Data;
                catch
                    disp('Cannot read trial c3d file')
                end
                
                fsampM = c3dStruct.Headers.VideoFrameRate; % camera frame rate
                AVRatio = c3dStruct.Headers.AVRatio; % camera to analogs ratio
                TotalFrames = c3dStruct.Headers.LastFrame; % total camera frames
                fsampK = AVRatio * fsampM; % forceplate sampling rate
                c3dTime = [0:length(c3d_Data.EMG)-1]./fsampK;
                otbTime = [0:length(otbEMG)-1]./fsamp;
                c3dEMG = c3d_Data.EMG;
                
                try
                    try
                        % synchronize forceplates and plot synchronized EMG
                        % signals
                        tic
                        [~, timeDiff, autoMatrix] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, fsampK, 0);
%                         aFig = figure(99);
%                         Autoax = axes('Parent',aFig);
%                         clf(Autoax)
%                         plot(Autoax, autoMatrix)
                        
                        TimeDelay = timeDiff;
                        zeroPadFront = zeros(timeDiff*2048,1);
                        
                        % upsample forceplates to fsamp
                        upsampPlates = resampler(cell2mat(struct2cell(c3d_Data.ForcePlates)'), fsampK, fsamp, 0);
                        TimeStop = (length(zeroPadFront) + size(upsampPlates,1)) / fsamp;
                        
                        [FPData, COPData] = bertec_COP(upsampPlates, 5, 2048);
                        % setup forceplate structure for easier data manipulation
                        fpLabels = fieldnames(c3d_Data.ForcePlates);
                        for ii = 1:size(fpLabels,1)
                            paddedFPVar = [zeroPadFront; upsampPlates(:,ii)]; % zero pads front of signal
                            paddedFPVar(length(paddedFPVar):trialLength)=0; % zero pads end of signal
                            ForcePlates.(fpLabels{ii}) = paddedFPVar;
                        end
                        
                        % upsample COP to fsamp
                        upsampCOP = resampler(cell2mat(struct2cell(c3d_Data.COP)'), fsampK, fsamp, 0);
                        
                        copLabels = {'LeftX', 'LeftY','RightX','RightY','WeightedX','WeightedY'};
                        %setup COP structure for easier data manipulation
                        copLabels = fieldnames(c3d_Data.COP);
                        for jj = 1:size(copLabels,1)
                            paddedCOPVar = [zeroPadFront; upsampCOP(:,jj)]; % zero pads front of signal
                            paddedCOPVar(length(paddedCOPVar) : trialLength) = 0; % zero pads end of signal
                            COP.(copLabels{jj}) = paddedCOPVar;
                        end
                        toc
                    catch
                        disp('no kinetics in this trials c3d file')
                    end
                    
                    try
                        tic
                        % synchronize Marker data
                        %[fsampMars, ~, ~] = emg_sync_v3(otbEMG, c3dEMG, fsamp, FrameRate, 1);
                        % upsample marker trajectories to fsamp
                        upsampMarkers = resampler(cell2mat(struct2cell(c3d_Data.Markers)'), fsampM, fsamp, 0);
                        
                        zeroFront = zeros(timeDiff * 2048, 1);
                        markerLabels = fieldnames(c3d_Data.Markers);
                        markerNumber = size(markerLabels,1);
                        zeroFrontMarkers = repmat(zeroFront,1, markerNumber*3);
                        FrontPadMarkers = [zeroFrontMarkers;upsampMarkers];
                        
                        padBackMarkers = zeros(length(size(FrontPadMarkers,1)+1:trialLength), 3);
                        
                        %setup Marker structure for easier data manipulation
                        loop = 0;
                        for kk = 1:3:(markerNumber*3)
                            loop = loop+1;
                            frontMarker = FrontPadMarkers(:, kk:kk+2);
                            fullMarkers = [frontMarker;padBackMarkers];
                            Markers.(markerLabels{loop}) = fullMarkers;%upsampMarkers(m1:m2, :); % zero pads front of signal
                        end
                        toc
                    catch
                        disp('no markers in this trials c3d file')
                    end
                    
                catch
                    % c3d data not processed
                end
                
                try
                    save(fullfile(kkFolder, save_KKFile), 'c3d_Data', 'Markers', 'fsampM', 'COP', 'ForcePlates', 'fsampK', 'c3dEMG', 'autoMatrix', 'AVRatio', 'TotalFrames', 'TimeDelay', 'TimeStop', '-v7.3','-nocompression');
                catch
                    save(fullfile(kkFolder, save_KKFile), 'c3d_Data', 'COP', 'ForcePlates','fsampK', 'c3dEMG','AVRatio', 'TotalFrames', 'TimeDelay', 'TimeStop', '-v7.3','-nocompression');
                end
            else
                % not processing c3d data if it already exists
            end
        end
        
    end
    
end

end
