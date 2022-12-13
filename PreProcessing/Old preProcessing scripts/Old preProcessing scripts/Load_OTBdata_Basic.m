
function [] = Load_OTBdata_Basic(fullFilePath, matrix, totalchan, matrixnames, trial, fsamp)

% DESCRIPTION
%
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
%
%   Implemented file structure for data processing; threshold analysis
%   moved to 'plotMC'; reduced memory overhead by eliminating 'fopen'
%   command; remove saturated channels automatically
%
% UPDATED   Matt Topley
% AUTHOR    Laura Miller
%           Christopher Thompson
%
% DATE CREATED      6-Feb-2013
%
% DATE MODIFIED   22-May-2013
% DATE MODIFIED   06-Sept-2022

%   To do: export more auxiliary channels; downsample force/auxiliary;
%   combine with 'SelectMCSeg' (or add this first);

[fileDir, trialname] = fileparts(fullFilePath);
clusterDir = fullfile(fileDir,'to_cluster');
trialFolder = fullfile(fileDir,trialname);
if exist(trialFolder, 'dir') == 7
else
    mkdir (trialFolder); % create unzip directory
end
fullTrialPath = fullfile(trialFolder,[trialname,'.otb']);
copyfile(fullFilePath,trialFolder);

unzip(fullTrialPath, trialFolder)
delete(fullTrialPath)

[allbadchan, data] = Preprocess_OTB_EMG(fullTrialPath,matrix,totalchan,fsamp,matrixnames);

if allbadchan == 'n'
    return
else
    allbadchanauto = 1:totalchan;
    allbadchanauto(matrix{1}) = [];
    allbadchan = [allbadchan; allbadchanauto'];
    
    EMGall = data(matrix{1},:);
    
    goodchan=[];
    for i = 1:(totalchan)
        if isempty(find(i==allbadchan))
            goodchan = [goodchan; i];
        end
    end
    
    if max(matrix{1}) == 64 && min(matrix{1}) == 1
        muscle = 'TA'
        
    elseif max(matrix{1}) == 128 && min(matrix{1}) == 65
        muscle = 'MG'
        
    elseif max(matrix{1}) == 192 && min(matrix{1}) == 129
        muscle = 'LG'
        
    elseif max(matrix{1}) == 256 && min(matrix{1}) == 193
        muscle = 'Sol'
        
    elseif max(matrix{1}) == 320 && min(matrix{1}) == 257
        muscle = 'dBF'
        
    elseif max(matrix{1}) == 384 && min(matrix{1}) == 321
        muscle = 'pBF'  
    else
        disp('Check your channels')
    end
    
    if exist(clusterDir,'dir') == 7
    else
        mkdir(clusterDir)
    end
    
    %%% save filename
    [~,subjectID] = fileparts(fileDir);
    saveid = strcat(subjectID,'_', trialname(17:21),'_',muscle,'_',num2str(trial), '.mat');
    
    % % %%% Quattro Setup Normal
    EMG = data(goodchan,:); % 256 no feedback
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
   
    save(fullfile(clusterDir,saveid),'EMG','EMGall','fsamp','allbadchan', 'Torque', 'Velocity', 'Position', 'Trigger', 'EMGTorque', ...
        'JR3Fx','JR3Fy','JR3Fz','JR3Mx','JR3My','JR3Mz', 'TargetEMG', 'TargetTorque', '-v7.3','-nocompression');
   
end
