
% function Swarm_Batch
clearvars('-except', 'f1Fig', 'ax1', 'ax2', 'ax3', 'ax4', 'f2ax', 'f2Fig') %%close all

keep = 0;
jointAction = 'DF';
feedback = 'Torque Feedback';
figFB = erase(feedback,' Feedback');
figIndex = 1;
figColor = 1;

fig2X = 11;
fig2Y = 12;

rootDir = cd;
filenames = dir(fullfile(rootDir, '**', '*TA_*.xls'));

sub2Remove = ~contains({filenames.name}, {'CoCoTest08','CoCotest010', 'CoCoTest_011'});
filenames(sub2Remove) = [];

folder2Remove = contains({filenames.folder}, {'\test\', 'Delta f'}) ;
filenames(folder2Remove) = [];

contraction2Remove = ~contains({filenames.folder}, feedback) ;
filenames(contraction2Remove,:) = [];

muReTime = [];
purse = [];
muDR = [];
for i = 1:size(filenames,1)
    try
        [RecruitmentTime, RecruitmentDR, percentMVC] = RecruitmentEMGamp(fullfile(filenames(i).folder,filenames(i).name));
        muReTime = [muReTime;RecruitmentTime];
        muDR = [muDR;RecruitmentDR];
        purse = [purse;percentMVC];
    end 
end

highUnits = muDR(purse>=5);
numbHigh = num2str(size(highUnits,1));
highUnitThresh = purse(purse>=5);
lowUnits = muDR(purse<5);
numbLow = num2str(size(lowUnits,1));
lowUnitThresh = purse(purse<5);

if exist('f2Fig') == 0
f2Fig = figure(2);%('Visible', 'off');
t1 = tiledlayout(2,2);
end 
f2Fig = figure(2)
nexttile(figIndex)
title(figFB)
legend('-DynamicLegend', 'Location', 'best');
hold on
%scatter(purse, muDR, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', figColor, 'DisplayName', figFB)
scatter(lowUnitThresh, lowUnits, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', figColor, 'DisplayName', [numbLow, ' Low Threshold Units'])
scatter(highUnitThresh, highUnits, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', figColor+1, 'DisplayName', [numbHigh, ' High Threshold Units'])
fitMod = polyfit(purse, muDR,1);
modFitline = polyval(fitMod,purse);
DisplayName = sprintf(['Slope ' num2str(fitMod(1))]); 
plot(purse,modFitline, 'SeriesIndex', figColor, 'LineWidth', 3,'DisplayName', DisplayName);
xlim([0 20]);
xlabel('Percentage of MVC EMG Amplitude');
ylim([0 35]);
ylabel('Recruitment DR');
hold off

if keep == 1
 t2Save = get(f2Fig,'Children');
    t2Save.Title.String = 'Recruitment threshold for % of EMG';
    set(f2Fig, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
    set(f2Fig, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
    saveas(f2Fig, fullfile('F:\Toolkit\Mirror\Delta f\DF Trials','TA_RecruitmentThreshold.pdf'))
end 

function [RecruitmentTime, RecruitmentDR, percentMVC] = RecruitmentEMGamp(fullFilename)

[folder, filename] = fileparts(fullFilename);
[fileDir,~] = fileparts(folder);
Under = strfind(filename,'_');

decompFile = [filename(1:Under(end-2)-1),'.mat'];
decompStruc = load(fullfile(fileDir,decompFile));
  
  if contains(filename,'CoCoTest_011')
      
      clusFile = [filename(1:(Under(5)-1)),'.mat'];
      load(fullfile(fileDir,clusFile));
      
      startRamp = str2num(filename(Under(7)+1:Under(8)-1));
      endRamp = str2num(filename(Under(8)+1:Under(9)-1));
  else
      clusFile = [filename(1:(Under(4)-1)),'.mat'];
      load(fullfile(fileDir,clusFile));
      
      startRamp = str2num(filename(Under(6)+1:Under(7)-1));
      endRamp = str2num(filename(Under(7)+1:Under(8)-1));
  end

Trialcut = [startRamp endRamp];
MUTime = [startRamp endRamp];
TQTime = [(startRamp*fsamp)+1 endRamp*fsamp];
timeTorque = TQTime(1)/fsamp:1/fsamp:MUTime(2);

MUPulses = decompStruc.MUPulses;
fsamp = decompStruc.fsamp;

MUFiring = SortUnits(MUPulses);
TRUEMU = [1:length(MUFiring)];

[trueMUFiring,trueMU] = TrimEdges2Ramp(MUTime,TRUEMU,MUFiring);


try
    TorqueFeedback = decompStruc.TorqueFeedback;
catch
    TorqueFeedback = decompStruc.JR3Mz;
end

try
    EMGFeedback = decompStruc.EMGFeedback;
catch
    EMGFeedback = decompStruc.EMGTorque;
end


EMGFeedback = abs(sum(EMGFeedback));
cutEMG = EMGFeedback(TQTime(1):TQTime(2));
newEMG = (cutEMG-(mean(cutEMG(1:fsamp))));
fs = 1/fsamp;
fn = fs/2;
[b,a] = butter(4,fn,'low');
y = filtfilt(b,a,newEMG);
realEMGFeedback = ((y - mean(y(1:100)))./(abs((max(y)-min(y)))))*10;

rampFig = figure(99)
plot(timeTorque,realEMGFeedback)
title(clusFile)
pause(1)


cutTorque = TorqueFeedback(TQTime(1):TQTime(2));
newTorque = (cutTorque-(mean(cutTorque(1:fsamp))));


if mean(newTorque)<0
    newTorqueFeedback = -newTorque./100;
else 
    newTorqueFeedback = newTorque./100;
end

RecruitmentTime =[];
RecruitmentDR = [];
percentMVC = [];

% figure(1)
% plot(timeTorque,realEMGFeedback)
for i = 1:length(trueMUFiring)
%     hold on
%     plot(trueMUFiring{i}(2:end),1./diff(trueMUFiring{i}(:)),'.')
%     plot(timeTorque(timeTorque(:) == trueMUFiring{i}(1)),realEMGFeedback(timeTorque(:) == trueMUFiring{i}(1)),'o')
    RecruitmentTime = timeTorque(timeTorque(:) == trueMUFiring{i}(1));
    RecruitmentDR = [RecruitmentDR;mean(1./diff(trueMUFiring{i}(1:3)))];
    percentMVC = [percentMVC;realEMGFeedback(timeTorque(:) == trueMUFiring{i}(1))];
end 
% hold off

% figure(1)
% plot(timeTorque,newTorqueFeedback)
% for i = 1:length(trueMUFiring)
%     hold on
%     plot(trueMUFiring{i}(2:end),1./diff(trueMUFiring{i}(:)),'.')
%     plot(timeTorque(timeTorque(:) == trueMUFiring{i}(2)),newTorqueFeedback(timeTorque(:) == trueMUFiring{i}(2)),'o')
%     RecruitmentTime = timeTorque(timeTorque(:) == trueMUFiring{i}(2));
%     RecruitmentDR = [RecruitmentDR;mean(1./diff(trueMUFiring{i}(1:3)))];
%     percentMVC = [percentMVC;newTorqueFeedback(timeTorque(:) == trueMUFiring{i}(2))];
% end 
% hold off
close(rampFig)

end


