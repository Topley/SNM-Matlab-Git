% function Swarm_Batch
clearvars('-except', 'f1Fig', 'ax1', 'ax2', 'ax3', 'ax4', 'f2ax', 'f2Fig') %%close all

rootdir = 'G:\Coco Subject Data\';
deltafFiles = dir(fullfile(rootdir, '**', '*Ramp*.pdf'));

[num,txt,raw] = xlsread('Master Subject Data Sheet.xls') ;

subject = 'Coco03';

findFiles = cellfun(@num2str,raw,'un',0);
[filenames] = getCocoFiles(findFiles, 'subject', subject, 'muscle', 'TA', 'feedback', 'Ramp', 'rampNum', 0, 'remove', 0, 'contraction', 'DF')

getCleanFiles = ~contains({deltafFiles.name},filenames(:,1));
deltafFiles(getCleanFiles,:) = [];


keeper1 = 0;
jointAction = 'DF';
feed = 'Coco';
feedback = [feed,' Feedback'];
rampNum = 3;
c = 5;
trialtype = feed;%['Ramp ', num2str(rampNum)];
fig2X = 24;
fig2Y = 24;

rootDir = cd;
% filenames = dir(fullfile(rootDir, '**', '*TA_*.xls'));
 
% contraction2Remove = ~contains({filenames.folder}, 'Delta f') ;
% filenames(contraction2Remove) = [];

joint2Remove = ~contains({filenames.folder}, feedback) ;
filenames(joint2Remove) = [];

contraction2Remove = contains({filenames.folder}, '\Test Subjects\') ;
filenames(contraction2Remove) = [];

sub2Remove = contains({filenames.name}, 'CoCoTest09') ;
filenames(sub2Remove) = [];

if exist('f1Fig') == 0 
f1Fig = figure(1);%('Visible', 'off');
    t = tiledlayout(2,2);
end 
f1Fig = figure(1)

allDFs = [];
allDRs = [];

keeper1 = 0;
allrawdfs = [];
for i = 1:size(filenames,1)
    try
        [rawdfs,avgdfs, data_MU] = rampSwarm(fullfile(filenames(i).folder,filenames(i).name));
        allDFs = [allDFs;avgdfs];
        allDRs = [allDRs; data_MU];
        allrawdfs = [allrawdfs; rawdfs];
    catch
    end
end


% allDFs(allDFs(:,9)>13,:) = [];
% allrawdfs(allrawdfs(:,11)>=4.5,:) = [];
% allDFs(allDFs(:,11)>=4.5,:) =[]; %& allDFs(:,11)<6,:) = [];
% allDRs(allDRs(:,11)>=4.5,:) = [];
% % % allDRs(allDRs(:,11)>=4 & allDRs(:,11)<6,:) = [];
%modulation = allDRs(:,12)-allDRs(:,11);
totalDF = size(allDFs, 1);
% binnumDF = round(sqrt(size(allDFs, 1)));
% [counts, edges] = histcounts(allDFs(:,3),binnumDF);
% binwidth = edges(2)-edges(1);
% bincenter = edges(1:end-1)+binwidth/2;
totalMU = size(allDRs, 1);
deltaf_binwidth = [-20:0.5:20];
reDR_binwidth = [0:0.2:10];
maxDR_binwidth = [0:0.5:30];
dereDR_binwidth = [0:0.2:10];
    
ax1 = nexttile(1)
hold on
legend('-DynamicLegend', 'Location', 'best');
dname = sprintf([trialtype ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDFs(:,3)); std(allDFs(:,3)); totalDF]);
s1 = swarmchart(ax1,rampNum*ones(1,length(allrawdfs(:,3))),allrawdfs(:,3), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c,'DisplayName', dname)
xlim([0 4]);
xticks([1 2 3 4]);
xticklabels({'Torque','EMG','Coco'});
ylim([-5 10]);
%xlabel('Ramp Number');
ylabel('Delta-F');

ax2 = nexttile(2)
hold on
legend('-DynamicLegend', 'Location', 'best');
dname = sprintf([trialtype ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDRs(:,11)); std(allDRs(:,11)); totalMU]);
s2 = swarmchart(ax2,rampNum*ones(1,length(allDRs(:,11))),allDRs(:,11),'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c,'DisplayName', dname)
xlim([0 4]);
xticks([1 2 3 4]);
xticklabels({'Torque','EMG','Coco'});
ylim([0 10]);
%xlabel('Ramp Number');
ylabel('Recruitment DR');

ax3 = nexttile(3)
hold on
legend('-DynamicLegend', 'Location', 'best');
dname = sprintf([trialtype ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDRs(:,12)); std(allDRs(:,12)); totalMU]);
s3 = swarmchart(ax3,rampNum*ones(1,length(allDRs(:,12))),allDRs(:,12), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c,'DisplayName', dname)
xlim([0 4]);
xticks([1 2 3 4]);
xticklabels({'Torque','EMG','Coco'});
ylim([5 25]);
%xlabel('Ramp Number');
ylabel('Max DR');

ax4 = nexttile(4)
hold on
legend('-DynamicLegend', 'Location', 'best');
dname = sprintf([trialtype ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDRs(:,13)); std(allDRs(:,13)); totalMU]);
Hysteresis = allDRs(:,11)-allDRs(:,13);
s4 = swarmchart(rampNum*ones(1,length(allDRs(:,13))),allDRs(:,13), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c,'DisplayName', dname)
%s4 = swarmchart(ax4, rampNum*ones(1,length(Hysteresis)),Hysteresis, 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c,'DisplayName', dname)
xlim([0 4]);
xticks([1 2 3 4]);
xticklabels({'Torque','EMG','Coco'});
ylim([0 10]);
%xlabel('Ramp Number');
ylabel('DeRecruitment DR');

if exist('f2Fig') == 0
f2Fig = figure(2);%('Visible', 'off');
t1 = tiledlayout(2,2);
end 
f2Fig = figure(2)
nexttile(rampNum)
title(feedback)
legend('-DynamicLegend', 'Location', 'best');
hold on
scatter(allDRs(:,fig2X), allDRs(:,fig2Y), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SeriesIndex', c, 'DisplayName', feed)
fitMod = polyfit(allDRs(:,fig2X),allDRs(:,fig2Y),1);
modFitline = polyval(fitMod,allDRs(:,fig2X));
DisplayName = sprintf(['Slope ' num2str(fitMod(1))]); 
plot(allDRs(:,fig2X),modFitline, 'SeriesIndex', c, 'LineWidth', 3,'DisplayName', DisplayName);
xlim([0 10]);
xlabel('Recruitment DR');
ylim([0 20]);
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

function [rawdfs, avgdfs, data_MU] = rampSwarm(fullFilename)

[folder, filename] = fileparts(fullFilename);
%[~,coherLoc] = fileparts(folder);

output = load(fullFilename);

outputclean = output;

outputclean(outputclean(:,1) == 0,:) = [];

outputclean(outputclean(:,4)<1,:) = [];
outputclean(outputclean(:,16)-outputclean(:,18)<1.5,:) = [];

%Replacing all unit comparisons with the average delta f for each unit
%rawdfs = outputclean;
rawdfs = [];
avgdfs = [];
Units = unique(outputclean(:,2));
for jk = 1:length(Units)
    placeholder = [];
    placeholder = outputclean(outputclean(:,2) == Units(jk),:);
    if size(placeholder,1) < 4
        %avgdfs(jk,:) = placeholder;
        avgdfs = [avgdfs;zeros(1,size(output,2))];
        rawdfs = [rawdfs;zeros(1,size(output,2))];
    else
        avgdfs(jk,:) = mean(placeholder);
        rawdfs = [rawdfs;placeholder];
    end
end

 avgdfs(avgdfs(:,1) == 0,:) = [];
  rawdfs(rawdfs(:,1) == 0,:) = [];
    data_MU = output;
    data_MU(data_MU(:,1) ~= 0,:) = [];
    data_MU(data_MU(:,2) == 0,:) = [];
end

