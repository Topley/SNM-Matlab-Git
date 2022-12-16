
function [] = Format_SNMLab_Data(MuscleList, TotalChans, AuxChans, varargin)

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
fsamp = 2048;
ArrayNumber = size(MuscleList,2);

%Parsing through files
for i = 1:nargin-3
    fileLoop = varargin{i};
    
    if contains(class(fileLoop), 'MatFile') || isempty(fileLoop)
        staticS = varargin{i};
        
    else
        [~, ~, ext] = fileparts(fileLoop);
        
        switch ext
            case '.otb'
                FileNamePathOTB = varargin{i};
            case '.otb+'
                FileNamePathOTB = varargin{i};
            case '.c3d'
                FileNamePathC3D = varargin{i};
        end
        
    end
end
%% Make trial directories if not created

% get OTB file information
[fileDir, otbTrialname] = fileparts(FileNamePathOTB);

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
    
    Muscle = MuscleList{ij}
    if isempty(Muscle)
        continue
    else
        LastChan = ij * 64; % calculate the last channel of the current array
        FirstChan = LastChan - 63; % calculate the first channel of the current array
        ArrayChannels = {FirstChan:LastChan}; % create cell array of the current array channels
        
        % preprocess OTB data
        [allbadchan, data] = Preprocess_OTB_EMG(ChannelMatrix, Muscle, ArrayChannels, TotalChans, fsamp);
        
        if ij == 1
            TrialNum = input('Enter the trial type (e.g., mvic = 0, ramp = 2): ', 's');
            if isempty(TrialNum) || length(TrialNum) < 1
                TrialNum = input('Missing trial type, Enter the trial type: ');%[60 85];
            end
        end
        
        % create to_cluster and kk filenames for saving
        [~,subjectID] = fileparts(fileDir);
        if ~contains(FileNamePathOTB, '_') && contains(FileNamePathOTB, '.otb+')
            newName = [otbTrialname,'00'];
            saveid = strcat(subjectID, '_', newName(end-4:end), '_', Muscle, '_', TrialNum, '.mat');
        elseif contains(otbTrialname, '_')
            newName = erase(otbTrialname,'_');
            saveid = strcat(subjectID, '_', newName(end-4:end), '_', Muscle, '_', TrialNum, '.mat');
        else
            saveid = strcat(subjectID, '_', otbTrialname(end-4:end), '_', Muscle, '_', TrialNum, '.mat');
        end
        
        % mvic = 0, quiet stance = 1, forward lean = 2, backward lean = 3, emg PF = 4, emg DF = 5, markerF = 6, markerB = 7, los = 8
        %% format EMG data for saving
        if allbadchan == 'n'
            continue
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
                Torque = ChannelMatrix(TotalChans-12,:);% 257 correct 3
                Velocity = ChannelMatrix(TotalChans-11,:);% 258 incorrect 4 60? possibly 180-120(JTAngle)
                Position = ChannelMatrix(TotalChans-10,:);% 259 incorrect 5 23960 ??
                Trigger = ChannelMatrix(TotalChans-9,:);% 260 correct 6
                EMGTorque = ChannelMatrix(TotalChans-8,:);% 261 EMG 7
                JR3Fx = ChannelMatrix(TotalChans-7,:);% 262 correct 8
                JR3Fy = ChannelMatrix(TotalChans-6,:);% 263 correct 9
                JR3Fz = ChannelMatrix(TotalChans-5,:);% 264 correct 10
                JR3Mx = ChannelMatrix(TotalChans-4,:);% 265 correct 11
                JR3My = ChannelMatrix(TotalChans-3,:);% 266 correct 12
                JR3Mz = ChannelMatrix(TotalChans-2,:);% 267 correct 13
                TargetEMG = ChannelMatrix(TotalChans-1,:);% 268 no feedback 14
                TargetTorque = ChannelMatrix(TotalChans,:); % 269 no feedback 15
                saveCluster = fullfile(toClusterFolder,saveid);
                
                try
                    save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'allbadchan',  'Torque', 'Velocity', 'Position', 'Trigger', 'EMGTorque', ...
                        'JR3Fx','JR3Fy','JR3Fz','JR3Mx','JR3My','JR3Mz', 'TargetEMG', 'TargetTorque', '-v7.3','-nocompression');
                catch
                    disp(['Unable to save Biodex setup for ', saveid])
                end
                
            elseif AuxChans < 13 && AuxChans > 1
                for iii = 0:AuxChans-1
                    Aux{iii} = ChannelMatrix(TotalChans-iii,:);
                end
                save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'Aux', 'allbadchan', '-v7.3', '-nocompression');
            else
                otbEMG = ChannelMatrix(TotalChans,:);
                save(fullfile(toClusterFolder,saveid),'EMG', 'EMGall', 'fsamp', 'otbEMG', 'allbadchan', '-v7.3', '-nocompression');
            end
        end
        
        %% read c3d file
        % if kk file does not already exist and there is c3d data available
        save_KKFile = replace(saveid, Muscle, 'kk'); % save kk file
        trialLength = length(otbEMG);
        
        if ij == 1 && isempty(FileNamePathC3D)
            continue
        elseif ij == 1 && ~isempty(FileNamePathC3D)
            try
                c3dStruct =  Read_C3D(FileNamePathC3D);
                c3d_Data = c3dStruct.Data;
                fsampM = c3dStruct.Headers.VideoFrameRate; % camera frame rate
                AVRatio = c3dStruct.Headers.AVRatio; % camera to analogs ratio
                TotalFrames = c3dStruct.Headers.LastFrame; % total camera frames
                fsampK = AVRatio * fsampM; % forceplate sampling rate
                c3dTime = [0:length(c3d_Data.EMG)-1]./fsampK;
                otbTime = [0:length(otbEMG)-1]./fsamp;
                c3dEMG = c3d_Data.EMG;
            catch
                disp('Cannot read trial c3d file')
            end
            
            try
                % synchronize forceplates
                tic
                [~, timeDiff, autoMatrix] = emg_sync_v3(otbEMG, c3d_Data.EMG, fsamp, fsampK, 1);
                fig2pos = get(figure(101),'position');
                fig2pos(1) = fig2pos(1)./2;
                set(figure(99),'position', fig2pos);
                TimeDelay = timeDiff;
                zeroPadFront = zeros(timeDiff*2048-1,1);
                upsampPlates = resampler(cell2mat(struct2cell(c3d_Data.ForcePlates)'), fsampK, fsamp, 0);
                TimeStop = (length(zeroPadFront) + size(upsampPlates,1)) / fsamp;
                fpLabels = fieldnames(c3d_Data.ForcePlates);
                reflect = size(zeroPadFront,1);
                for ii = 1:size(fpLabels,1)
                    reflectPad = flip(upsampPlates(1:reflect,ii));
                    paddedFPVar = [reflectPad; upsampPlates(:,ii)] - mean(reflectPad(1:100)); % zero pads front of signal
                    ForcePlates.(fpLabels{ii}) = paddedFPVar(1:trialLength);
                end
                
                upsampCOP = resampler(cell2mat(struct2cell(c3d_Data.COP)'), fsampK, fsamp, 0);
                copLabels = {'LeftX', 'LeftY','RightX','RightY','WeightedX','WeightedY'};
                copLabels = fieldnames(c3d_Data.COP);
                for jj = 1:size(copLabels,1)
                    reflectPad = flip(upsampCOP(1:reflect,jj));
                    paddedCOPVar = [reflectPad; upsampCOP(:,jj)]- mean(reflectPad(1:100)); % zero pads front of signal
                    COP.(copLabels{jj}) = paddedCOPVar(1:trialLength);
                end
                toc
            catch
                disp('no kinetics in this trials c3d file')
                ForcePlates = [];
                COP = [];
            end
            
            try
                tic
                % synchronize Marker data
                upsampMarkers = resampler(cell2mat(struct2cell(c3d_Data.Markers)'), fsampM, fsamp, 0);
                zeroFront = timeDiff * 2048 - 1;
                markerLabels = fieldnames(c3d_Data.Markers);
                markerNumber = size(markerLabels,1);
                reflectMar = upsampMarkers([1:zeroFront], :);
                FrontPadMarkers = [flip(reflectMar);upsampMarkers];
                
                loop = 0;
                for kk = 1:3:(markerNumber*3)
                    loop = loop+1;
                    fullMarkers = FrontPadMarkers(:, kk:kk+2);
                    %fullMarkers = [frontMarker;padBackMarkers];
                    Markers.(markerLabels{loop}) = fullMarkers(1:trialLength, :);%upsampMarkers(m1:m2, :); % zero pads front of signal
                end
                toc
                
                try
                    [JointAngle, JointCenter] = ProcessKinematics('right', staticS, Markers);
                end
                
                try
                    [JointAngle, JointCenter] = ProcessKinematics('left', staticS, Markers, JointAngle, JointCenter);
                end
                
            catch
                disp('no markers in this trials c3d file')
                Markers = [];
                
                if ~exist('JointAngle', 'var')
                    JointAngle = [];
                end
                
                if ~exist('JointCenter', 'var')
                    JointCenter = [];
                end
            end
            
            save(fullfile(kkFolder, save_KKFile), 'c3d_Data', 'Markers', 'JointAngle', 'JointCenter', 'COP', 'ForcePlates', 'c3dEMG', 'autoMatrix', 'TimeDelay', 'TimeStop', '-v7.3','-nocompression');
        else
            
        end
    end
end
end
