% function Swarm_Batch
clearvars('-except', 'f1Fig', 'ax1', 'ax2', 'ax3', 'ax4', 'f2ax', 'f2Fig') %%close all

keeper1 = 0;
jointAction = 'DF';
feed = 'Torque';
feedback = [feed,' Feedback'];
rampNum = 3;
c = 5;
trialtype = feed;%['Ramp ', num2str(rampNum)];
fig2X = 24;
fig2Y = 24;

rootDir = cd;
filenames = dir(fullfile(rootDir, '**', '*TA_*.xls'));
 
% contraction2Remove = ~contains({filenames.folder}, 'Delta f') ;
% filenames(contraction2Remove) = [];

joint2Remove = ~contains({filenames.folder}, feedback) ;
filenames(joint2Remove) = [];

% contraction2Remove = ~contains({filenames.folder}, '\test\') ;
% filenames(contraction2Remove) = [];

sub2Remove = contains({filenames.name}, 'CoCoTest09') ;
filenames(sub2Remove) = [];

for i = 1:size(filenames,1)
    try
        [DR,EMGamp] = getEMGvsDR(fullfile(filenames(i).folder,filenames(i).name));
      
    catch
    end
end

% if exist('f2Fig') == 0
% f2Fig = figure(2);%('Visible', 'off');
% t1 = tiledlayout(2,2);
% end 
% f2Fig = figure(2)
% nexttile(rampNum)
% title(feedback)
% legend('-DynamicLegend', 'Location', 'best');
% hold on
% scatter(allDRs(:,fig2X), allDRs(:,fig2Y), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', c, 'DisplayName', feed)
% fitMod = polyfit(allDRs(:,fig2X),allDRs(:,fig2Y),1);
% modFitline = polyval(fitMod,allDRs(:,fig2X));
% DisplayName = sprintf(['Slope ' num2str(fitMod(1))]); 
% plot(allDRs(:,fig2X),modFitline, 'SeriesIndex', c, 'LineWidth', 3,'DisplayName', DisplayName);
% xlim([0 10]);
% xlabel('Recruitment DR');
% ylim([0 20]);
% ylabel('Max DR');
% hold off

function [DR,EMGamp] =  getEMGvsDR(fullFilename)
fsamp =2048;
[folder, filename] = fileparts(fullFilename);
[fileDir,~] = fileparts(folder);
Under = strfind(filename,'_');  
decompFile = [filename(1:Under(end-2)-1),'.mat'];
decompStruc = load(fullfile(fileDir,decompFile));

m1Pulses = decompStruc.MUPulses;
%%% sort units
m1MotorUnits = SortUnits(m1Pulses);

EMG = decompStruc.EMG;
absRawEMG = sum(abs(EMG));
    rectifiedEMG = absRawEMG-mean(absRawEMG(1:fsamp));
    rmsEMG = rms(rectifiedEMG, 500);
    rmsEMG = rmsEMG-mean(rmsEMG(1:100));
    EMGFeedback = movmean(rmsEMG, fsamp);
    fs = 1/fsamp;
    fn = fs/4
    wn = [1 100]
    [b,a] = butter(3,fs,'low')
    y = filtfilt(b,a,rmsEMG);
    rEMGFeedback = ((y - mean(y(1:100)))./(abs((max(y)-min(y)))))*10;

    plot([0:length(rEMGFeedback)-1]/fsamp,rEMGFeedback)
    1./diff(m1MotorUnits{1}(2:3))*fsamp
    for i =1:length(m1MotorUnits)
    hold on
        %plot(rEMGFeedback(m1MotorUnits{i}(1)),mean((1./diff(m1MotorUnits{i}(1:3)))*fsamp),'.', 'MarkerSize',20) 
        plot(m1MotorUnits{i}(1)/fsamp,mean((1./diff(m1MotorUnits{i}(1:3)))*fsamp),'.', 'MarkerSize',20)
    end



%Replacing all unit comparisons with the average delta f for each unit
%rawdfs = outputclean;
r
end 