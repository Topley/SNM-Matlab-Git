
close all;
clear all;

% Open Decomposed MU File
filename = uigetfile;
load(filename);
Tableaumap = Tableau;   % colormap if user doesn't have predefined colormap

Under = strfind(filename,'_');                      % search for file delimiters
clusFile = [filename(1:(Under(end-1)-1)),'.mat'];   % Find cluster file if neccessary

if exist(clusFile, 'file') == 2
    try
        load(clusFile, 'fsamp', 'Torque', 'JR3Mz');
    end
end

%%% Probably should use this and not RMSE
realTorque = ((JR3Mz - mean(JR3Mz(1:fsamp))));
newTorque = realTorque./100;
% plot([1:length(newTorque)-1]/fsamp,diff(newTorque*1000)
% hold on
% plot(MUPulses{end}(2:end)/fsamp,1./diff(MUPulses{end}).*fsamp,'o')

MUFiring = MUPulses;

firing =[];
for i = 1:length(MUFiring)
    MuF = MUFiring{i};
    index = round(MuF);
    firing(i,index) = 1;
end
firing(:,length(firing)+fsamp)=0; %adds a second of 0s to the end
AVERAGE = 1./sum(hanning(fsamp*2))*filtfilt(hanning(fsamp*2),1,firing')';%AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
[pks, locs] = findpeaks(mean(AVERAGE), 'NPeaks', 3,'MinPeakDistance',fsamp*15, 'MinPeakHeight', 2);
locs = locs./fsamp;
avgvector = mean(AVERAGE);
figure(1);clf
    hold on
    plot([0:length(newTorque)-1]/fsamp, newTorque)
    plot([0:length(AVERAGE)-1]/fsamp,mean(AVERAGE))
    plot(locs, pks,'*')
%     findpeaks([0:length(AVERAGE)-1]/fsamp,mean(AVERAGE), 'MinPeakDistance',20);
    hold off
  startramp = locs-10;
  stopramp = locs+10;
% plot(AVERAGE(end,startramp(1)*fsamp:stopramp(1)*fsamp))

% MUFiring = MUPulses;
% MUFiring(find(cellfun(@isempty,MUFiring))) = [];
% [~,FiringLen] = sort(cellfun(@length,MUFiring));
% MUFiring = MUFiring(FiringLen);


%%% start and end f Ramp
% Trialcut1 = input('contraction 1: '); %[5 30];
% Trialcut2 = input('contraction 2: '); %[35 60];
% Trialcut3 = input('contraction 3: '); %[65 90];

Trialcut1 = [startramp(1) stopramp(1)]; %[5 30];
Trialcut2 = [startramp(2) stopramp(2)]; %[35 60];
Trialcut3 = [startramp(3) stopramp(3)]; %[65 90];

MUTime1 = [Trialcut1(1) Trialcut1(2)];
TQTime1 = [(Trialcut1(1)*fsamp)+1 (Trialcut1(2)*fsamp)];

newTorqueFeedback1 = newTorque(TQTime1(1):TQTime1(2));
timeTorque1 = TQTime1(1)/fsamp:1/fsamp:MUTime1(2);

%second ramp
MUTime2 = [Trialcut2(1) Trialcut2(2)];
TQTime2 = [(Trialcut2(1)*fsamp)+1 (Trialcut2(2)*fsamp)];

newTorqueFeedback2 = newTorque(TQTime2(1):TQTime2(2));
timeTorque2 = TQTime2(1)/fsamp:1/fsamp:MUTime2(2);

% Third ramp
MUTime3 = [Trialcut3(1) Trialcut3(2)];
TQTime3 = [(Trialcut3(1)*fsamp)+1 (Trialcut3(2)*fsamp)];

newTorqueFeedback3 = newTorque(TQTime3(1):TQTime3(2));
timeTorque3 = TQTime3(1)/fsamp:1/fsamp:MUTime3(2);


RMSETF1 = sqrt((0-newTorqueFeedback1).^2);
normRMSETF1 = mean(RMSETF1)./(max(RMSETF1)-min(RMSETF1));

RMSETF2 = sqrt((0-newTorqueFeedback2).^2);
normRMSETF2 = mean(RMSETF2)./(max(RMSETF2)-min(RMSETF2));

RMSETF3 = sqrt((0-newTorqueFeedback3).^2);
normRMSETF3 = mean(RMSETF3)./(max(RMSETF3)-min(RMSETF3));

figure;
hold all

plot(timeTorque1, newTorqueFeedback1,'SeriesIndex', 2, 'DisplayName', ['RMSEfb ', num2str(normRMSETF1)])
plot(timeTorque2, newTorqueFeedback2,'SeriesIndex', 4, 'DisplayName', ['RMSEfb ', num2str(normRMSETF2)])
plot(timeTorque3, newTorqueFeedback3,'SeriesIndex', 6, 'DisplayName', ['RMSEfb ', num2str(normRMSETF3)])

titletxt1 = sprintf('Ramp 1 start %0.0f & stop %0.0f',startramp(1), stopramp(1));
titletxt2 = sprintf('Ramp 2 start %0.0f & stop %0.0f',startramp(2), stopramp(2));
titletxt3 = sprintf('Ramp 3 start %0.0f & stop %0.0f',startramp(3), stopramp(3));
title({titletxt1,titletxt2,titletxt3})

legend('show')

%sqrt(mean((0-newTorque1).^2))./(max(newTorque1)-min(newTorque1))