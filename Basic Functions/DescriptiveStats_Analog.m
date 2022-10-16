clearvars % better than clear all
%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls');
% subject = 'Coco01';
% subject = 'Coco03_23663_TA_5_v23decomposed';
[filenames] = getCocoFiles(masterSheet, 'muscle', 'TA', 'feedback', 'Ramp', 'contraction', 'Coco')

%%%% User needs to adjust these before running script
if contains(filename, 'Coco01')
    usermaxTorque = 52;
elseif contains(filename, 'Coco02')
    usermaxTorque = 53.5;
elseif contains(filename, 'Coco03')
    usermaxTorque = 45;
elseif contains(filename, 'Coco04')
    usermaxTorque = 46;
end

Auto = 1;                   % 1; %no bad MUs 1 or 0
torqueRMSE = [];
peakEMG = [];
peakTorq = [];
loop = 0;
for i = 1:size(filenames,1)
    subject = filenames{i}(1:6)
    for ij = 1:3
        loop = loop+1;
        warning('off', 'MATLAB:LOAD:VariableNotFound'); % disable warning if variable can't be loaded
        try
            [maxTorque, maxEMG, normRMSETorque] =  rampDescriptive(fullfile(rootDir,subject,'decomposed',filenames(i)),usermaxTorque, Auto,ij);
            peakTorq = [peakTorq,maxTorque];
            peakEMG = [peakEMG,maxEMG];
            torqueRMSE = [torqueRMSE,normRMSETorque];
        catch
            disp('Could not complete delta f for previous file');
            continue
        end
    end
end


function [maxTorque, maxEMG, normTorqueRMSE] = rampDescriptive(fullFilename,usermaxTorque, Auto, ramp)

close all;

%%% Set up directories - change if > 2 subdirectories
%%% Set up directories - change if > 2 subdirectories
[fileDir,filename,clusterFile, pdfdir] = setupDirectories(fullFileName, 'ramp');

% Open Decomposed MU File
load(fullfile(Feedbackdir,filename),'MUPulses', 'TraceFeedback')
%load(filename, 'MUPulses', 'Torque', 'TorqueFeedback', 'EMGFeedback');
Tableaumap = Tableau;   % colormap if user doesn't have predefined colormap

try
    load(fullfile(fileDir,clusterFile), 'fsamp', 'Torque');
    temp = load(fullfile(fileDir,clusterFile), 'JR3Mz');
    TorqueFeedback = temp.JR3Mz;
end

try
    EMGFeedback = [];
    load(fullfile(fileDir,clusterFile), 'EMG');
end

absRawEMG = abs(EMG(27,:));
rectifiedEMG = absRawEMG-mean(absRawEMG(1:fsamp));
rmsEMG = rms(rectifiedEMG, 500);
EMGFeedback = movmean(rmsEMG, fsamp);
EMGFeedback = EMGFeedback-mean(EMGFeedback(1:fsamp));
fs = 1/fsamp;
[b,a] = butter(3,fs,'low');
EMGthreshold = filtfilt(b,a,rmsEMG);

%%%%% load agonist units
if exist('fsamp') == 0
    fsamp = 2048;
end

checkTorque = Torque-mean(Torque(1:fsamp));
if abs(mean(checkTorque)./usermaxTorque) < 1
    coco = 1;
else
    coco = 0;
end


% Find antagonist file, adjust these if analyzing different muscles
if coco == 1;
    [AntagEMGFeedback] = findAntagonistEMG(fullfile(fileDir,clusterFile), 'TA', 'Sol');
end

%%%%%
[rampStart] = findRamps(TraceFeedback,EMGFeedback, coco);

if ramp == 1
    Trialcut = [rampStart(1) rampStart(1)+30];
elseif ramp == 2
    Trialcut = [rampStart(2) rampStart(2)+30];
else
    Trialcut = [rampStart(3) rampStart(3)+30];
end

if Trialcut(1) < 1
    Trialcut(1) = 1;
end
stax = Trialcut(1);

%%%% plot actual torque, torque feedback, and EMG feedback

%%%% cut trialcut length is ok
if Trialcut(2)>(length(EMGFeedback)/fsamp)-1
    Trialcut(2) = (floor(length(EMGFeedback)/fsamp))-1;
end

%%%%
%This is where MUFiring, torque, et al should be cut
MUTime = [Trialcut(1) Trialcut(2)];
try
    %%% cut real torque and Torque Feedback
    [TQTime, timeTorque, cutTorque, cutTorqueFeedback, pfTQfeedback, dfTQfeedback] = TrimTorque(MUTime,Torque, TorqueFeedback, usermaxTorque, coco);
    
    %%% cut EMG feedback
    EMGFeedback = ((EMGFeedback - mean(EMGFeedback(1:fsamp)))./(abs((max(EMGFeedback)-min(EMGFeedback)))))*20;
    %      EMGFeedback = ((rectifiedEMG - mean(rectifiedEMG(1:fsamp)))./(abs((max(rectifiedEMG)-min(rectifiedEMG)))))*20;
    newEMGFeedback = EMGFeedback(TQTime(1):TQTime(2));
    %     EMGthreshold = EMGthreshold(TQTime(1):TQTime(2));
    %     EMGthreshold = EMGthreshold -mean(EMGthreshold (1:fsamp));
    %     newEMGFeedback = rectifiedEMG(TQTime(1):TQTime(2));
    %%% cut antagnoist EMG feedback if applicable
    try
        AntagEMGFeedback = ((AntagEMGFeedback - mean(AntagEMGFeedback(1:fsamp)))./(abs((max(AntagEMGFeedback)-min(AntagEMGFeedback)))))*10;
        newAntagEMGFeedback = AntagEMGFeedback(TQTime(1):TQTime(2));
    end
catch
end

maxTorque = max(abs((cutTorqueFeedback./usermaxTorque).*100));
maxEMG = max(newEMGFeedback);
torqueRMSE = sqrt((cutTrace-cutTorqueFeedback).^2);
normTorqueRMSE = mean(torqueRMSE)./(max(torqueRMSE)-min(torqueRMSE));

end