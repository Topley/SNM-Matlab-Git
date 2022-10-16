% function Swarm_Batch
clearvars('-except', 'f1Fig', 'ax1', 'ax2', 'ax3', 'ax4', 'f2ax', 'f2Fig') %%close all

keeper1 = 0;
jointAction = 'DF';
feed = 'EMG';
feedback = [feed,' Feedback'];
rampNum = 2;
c = 3;
trialtype = feed;%['Ramp ', num2str(rampNum)];
fig2X = 11;
fig2Y = 12;

rootDir = cd;
filenames = dir(fullfile(rootDir, '**', '*pooled_coher*.mat'));

joint2Remove = ~contains({filenames.name}, 'TA') ;
filenames(joint2Remove) = [];

% sub2Remove = contains({filenames.name}, {'CoCoTest09', 'CoCotest03'}) ;
% filenames(sub2Remove) = [];

contraction2Remove = ~contains({filenames.folder}, feedback) ;
filenames(contraction2Remove) = [];

if exist('f1Fig') == 0 
f1Fig = figure(1);%('Visible', 'off');
    t = tiledlayout(2,2);
end 

maxDRs =[];
reDRs = [];
deReDRs = [];
for i = 1:size(filenames,1)
    try
        [maxDR, reDR, deReDR] = holdSwarm(fullfile(filenames(i).folder,filenames(i).name));
        maxDRs = [maxDRs;maxDR];
        reDRs = [reDRs;reDR];
        deReDRs = [deReDRs;deReDR];
    catch
    end
end
reDRs(maxDRs>35) = [];
deReDRs(maxDRs>35) = [];
maxDRs(maxDRs>35) = [];

if exist('f2Fig') == 0
f2Fig = figure(2);%('Visible', 'off');
t1 = tiledlayout(2,2);
end 
f2Fig = figure(2)
nexttile(rampNum)
title(feedback)
legend('-DynamicLegend', 'Location', 'best');
hold on
scatter(reDRs, maxDRs, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', c, 'DisplayName', feed)
fitMod = polyfit(reDRs,maxDRs,1);
modFitline = polyval(fitMod,reDRs);
DisplayName = sprintf(['Slope ' num2str(fitMod(1))]); 
plot(reDRs,modFitline, 'SeriesIndex', c, 'LineWidth', 3,'DisplayName', DisplayName);
xlim([0 20]);
xlabel('Recruitment DR');
ylim([0 35]);
ylabel('Max DR');
hold off

if keeper1 == 1
    tSave = get(f1Fig, 'Children');
    tSave.Title.String = 'All Collected Ramps';
    set(f1Fig, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
    set(f1Fig, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
    saveas(f1Fig, fullfile('F:\Toolkit\Mirror\Delta f\DF Trials','TA_Deltaf_Plots_4compsMin.pdf'))
    
    
    t2Save = get(f2Fig,'Children');
    t2Save.Title.String = 'Recruitment and rate modulation slopes';
    set(f2Fig, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
    set(f2Fig, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
    saveas(f2Fig, fullfile('F:\Toolkit\Mirror\Delta f\DF Trials','TA_RateModulationSlopes_4compsMin.pdf'))
end
% end

function [maxMU,reMU, deReMU] = holdSwarm(fullFilename)
fsamp = 2048;
[folder, filename] = fileparts(fullFilename);
Under = strfind(filename,'_');    
%[~,coherLoc] = fileparts(folder);
decompFile = [filename(1:Under(end-3)-1),'.mat'];
temp1 = load(fullfile(folder,decompFile), 'MUPulses');
m1Pulses = temp1.MUPulses;
reMU = [];
maxMU =[];
deReMU = [];
%%% sort units
m1MotorUnits = SortUnits(m1Pulses);
for i = 1:length(m1MotorUnits)
instantDR = 1./(diff(m1MotorUnits{i})).*fsamp;
reMU = [reMU;mean(instantDR(1:3))];
maxMU = [maxMU;mean(maxk(instantDR,3))];
deReMU = [deReMU;mean(instantDR(end-3:end))];
% plot(m1MotorUnits{i}(2:end)./fsamp, instantDR, 'o')
end 
end

