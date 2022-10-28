clear all;
close all;
clc;
%%%NORMAL%%%
[FileName,PathName] = uigetfile({'*.mat'},'Select the MATLAB Decomposed file');
load(fullfile(PathName,FileName));
ind = strfind(FileName,'_');

clusterPath = replace(PathName, 'decomposed', 'to_cluster');
kkPath = replace(PathName, 'decomposed', 'kk files');

kkFile = [FileName(1:ind(2)), 'kk',FileName(ind(3):ind(4)-1) '.mat'];
kkFile = replace(kkFile, 'kk_5', 'kk_6');
vars = matfile(fullfile(clusterPath,FileName(1:ind(end)-1)));
varCOP = load(fullfile(kkPath,kkFile));

cop = varCOP.COP.WeightedY;
unpadded = cop ~= 0;
padded = cop == 0;
offset = min(cop(unpadded));
cop(padded) = offset;

EMGall = vars.EMGall;
EMG = vars.EMG;

cop = ((cop - mean(cop(1:100)))./(abs((max(cop)-min(cop)))))*10;
% [~,EMGFeedback] = EMGpro(vars.EMGall, 'channel', 27, 'filter', {'filtRMS', 500, 250});
% EMGFeedback = EMGFeedback./10;


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

save([FileName(1:end-4),'_MUCLEANED.mat'],'-v7.3');

%%%NORMAL%%%
improve_decomposition_V15([PathName,FileName(1:end-4),'_MUCLEANED.mat'],'MUPulses','IPTs','SIL','EMG','SIG','FACTOR','PEREIG','fsamp','FREQ', 'cop',[0 length(EMGall)/2048],'LSD');

%load([FileName(1:end-4),'_MUCLEANED.mat']);

% MUPulses(cellfun(@isempty,MUPulses))=[];
% SIL(cellfun(@isempty,SIL))=[];
% IPTs(cellfun(@isempty,IPTs))=[];

clear trial PathName ans ROW COL ind FIRSTTIME vars
save([FileName(1:end-4),'_MUCLEANED.mat'],'-v7.3','-nocompression');

