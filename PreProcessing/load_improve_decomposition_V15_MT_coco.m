clear all;
close all;
clc;
%

% if length(dir('*XXX*'))>0
% disp('Use the program "MUCk2020_namechanger_RUNFIRST.m" to enter your initials')
% else

%%%NORMAL%%%
[FileName,PathName] = uigetfile({'*.mat'},'Select the MATLAB Decomposed file');
load(fullfile(PathName,FileName));
ind = strfind(FileName,'_');
try
    vars = matfile(fullfile(PathName,FileName(1:ind(end)-1)));
    EMG = vars.EMG; 
    Torque = vars.AnalogData(10,:);
catch
%    load(fullfile(PathName,FileName(1:ind(end)-1)),'EMG','TraceFeedback', 'Feedback','EMGTorque');
end


try
    TorqueFeedback = (JR3Mz - mean(JR3Mz(1:100)))./43;
catch
    TorqueFeedback = ((Torque - mean(Torque(1:100)))./(abs((max(Torque)-min(Torque)))))*10;
end

if mean(TorqueFeedback)<0
    TorqueFeedback = TorqueFeedback.*-1;
end
 
try
      rawEMG = abs(sum(EMG));
   
catch
  
    rawEMG = abs(EMGTorque-mean(EMGTorque(1:fsamp)));
end

rmsEMG = rms(rawEMG, 500);
EMGFeedback = movmedian(rmsEMG, fsamp);


pchan(1,:)=[nan 1 2 3 4];
pchan(2,:)=[16 8 7 6 5];
pchan(3,:)=[15 14 13 12 11];
pchan(4,:)=[19 18 17 9 10];
pchan(5,:)=[20 21 22 23 24];
pchan(6,:)=[28 29 30 31 32];
pchan(7,:)=[25 26 27 33 34];
pchan(8,:)=[35 36 37 38 39];
pchan(9,:)=[45 46 47 48 40];
pchan(10,:)=[44 43 42 41 49];
pchan(11,:)=[54 53 52 51 50];
pchan(12,:)=[55 56 64 63 62];
pchan(13,:)=[57 58 59 60 61];

try
    for ROW = 1:13
        for COL = 1:5
            if isnan(pchan(ROW,COL))
                SIG{ROW,COL} = [];
            else
                SIG{ROW,COL} = EMGall(pchan(ROW,COL),:);
            end
        end
    end
catch
    
    for ROW = 1:13;
        for COL = 1:5;
            SIG{ROW,COL} = zeros(1,length(EMG));
        end
    end
end
clear rmsEMG rawEMG 
save([FileName(1:end-4),'_MUCLEANED.mat'],'-v7.3');

%%%NORMAL%%%
improve_decomposition_V15([PathName,FileName(1:end-4),'_MUCLEANED.mat'],'MUPulses','IPTs','SIL','EMG','SIG','FACTOR','PEREIG','fsamp','FREQ','TorqueFeedback',[0 length(EMG)/2048],'LSD');

load([FileName(1:end-4),'_MUCLEANED.mat']);

% MUPulses(cellfun(@isempty,MUPulses))=[];
% SIL(cellfun(@isempty,SIL))=[];
% IPTs(cellfun(@isempty,IPTs))=[];

clear trial PathName ans ROW COL ind FIRSTTIME rmsEMG rmsEMG27 rawEMGFeedback rawEMG rawEMG27
save([FileName(1:end-4),'_MUCLEANED.mat'],'-v7.3','-nocompression');


