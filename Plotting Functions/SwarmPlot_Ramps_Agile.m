 Swarm_Batch pairwise and unitwise
clearvars('-except', 'f1Fig', 'ax1', 'ax2', 'ax3', 'ax4', 'f2ax', 'f2Fig') %%close all

%%%% Figure setup %%%
keeper1 = 0;
rampNum = 3;

if rampNum == 1
    feedback = 'DF Torque';
    c = 1;
    contract = 'Coco';
    trial = 'Ramp';
elseif rampNum == 2
    feedback = 'Variable Coco';
    c = 3;
    contract = 'Coco';
    trial = 'Ramp';
elseif rampNum == 3
    feedback = 'Invariable Coco';
    c = 5;
    contract = 'Coco';
    trial = 'Ramp';
end

%%%% File Lookup %%%%
rootDir = 'G:\Coco Subject Data\';
deltafFiles = dir(fullfile(rootDir, '**', '*Ramp3*.pdf'));

%%%% Pick files directory - use wildcard t select specific muscles
masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls'); 
% subject = 'Coco03_23663_TA_5_v23decomposed';
[filenames] = getCocoFiles(masterSheet,  'muscle', 'TA', 'feedback', trial, 'contraction', contract)

getCleanFiles = ~contains({deltafFiles.name},filenames(:,1));
deltafFiles(getCleanFiles,:) = [];
% scatter(allrawdfs(:,25),allrawdfs(:,3), 'filled', 'SeriesIndex', 1) 

% ramp2Remove = ~contains({deltafFiles.name}, ['Coco01_24941']);%_TA_5_v23decomposed_MUCLEANED_62_82_Ramp3']) ;
% deltafFiles(ramp2Remove) = [];

if exist('f1Fig') == 0
    f1Fig = figure(1);%('Visible', 'off');
    t = tiledlayout(2,2);
end
f1Fig = figure(1);

allDFs = [];
allDRs = [];
allrawdfs = [];
trial = [];
subNum = [];
for i = 1:size(deltafFiles,1)
%for i = 1:size(filenames,1)
    try
        subTrial = str2num(deltafFiles(i).name(6));
        conNum = str2num(deltafFiles(i).name(end-4));
        [rawdfs,avgdfs, data_MU] = rampSwarm1020(fullfile([deltafFiles(i).folder, '\results_excel\'],[deltafFiles(i).name(1:end-9),'deltaf.xls']));
        %[rawdfs,avgdfs, data_MU] = rampSwarm1020(fullfile(filenames(i).folder, filenames(i).name));
        allDFs = [allDFs;avgdfs];
        allDRs = [allDRs; data_MU];
        allrawdfs = [allrawdfs; rawdfs];
        numUnitwise = zeros(1,size(avgdfs,1));
        subUnitwise = zeros(1,size(avgdfs,1));
        numUnitwise(:) = conNum;
        subUnitwise(:) = subTrial;
        trial = [trial;numUnitwise'];
        subNum = [subNum;subUnitwise'];
    catch
        disp(['Could not process ', deltafFiles(i).name])
    end
end

findMin = [allDRs(:,11) allDRs(:,13)];
M = min(findMin,[],2);
totalDF = size(allDFs, 1);
totalRAW = size(allrawdfs, 1);
totalMU = size(allDRs, 1);
deltaf_binwidth = [-20:0.5:20];
reDR_binwidth = [0:0.2:10];
maxDR_binwidth = [0:0.5:30];
dereDR_binwidth = [0:0.2:10];

ax1 = nexttile(1);
hold on
s1 = swarmchart(ax1,rampNum*ones(1,length(allrawdfs(:,3))),allrawdfs(:,3), 'filled', 'MarkerFaceAlpha', 0.8,'MarkerFaceColor', [.8 .8 .8]);
% title('Pairwise Delta f');
% dname = sprintf([feedback ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allrawdfs(:,3)); std(allrawdfs(:,3)); totalRAW]);
% holdSubtitle = string(ax1.Subtitle.String);
% if holdSubtitle == ""
%     set(ax1.Subtitle,'String', dname)
% else
%     set(ax1.Subtitle,'String', [holdSubtitle;dname])
% end 
set(ax1,'TickDir','out');
xlim([0 4]);
xticks([]);
ax1.XColor = 'w';
ylim([-2 10]);
dname2 = sprintf([feedback ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDFs(:,3)); std(allDFs(:,3)); totalDF]);
s2 = swarmchart(ax1,rampNum*ones(1,length(allDFs(:,3))),allDFs(:,3),'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c, 'SizeData', 50);
title('Unitwise Delta f');
holdSubtitle2 = string(ax1.Subtitle.String);
if holdSubtitle2 == ""
    set(ax1.Subtitle,'String', dname2)
else
    set(ax1.Subtitle,'String', [holdSubtitle2;dname2])
end 

ax2 = nexttile(2);
hold on
dname = sprintf([feedback ',' ' %0.2f ' char(177) ' %0.2f ' '(%d Total)'], [mean(allDRs(:,12)); std(allDRs(:,12)); totalMU]);
s2 = swarmchart(ax2,rampNum*ones(1,length(allDRs(:,12))),allDRs(:,12), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c);
title('Max Discharge Rate');
holdSubtitle2 = string(ax2.Subtitle.String);
if holdSubtitle2 == ""
    set(ax2.Subtitle,'String', dname)
else
    set(ax2.Subtitle,'String', [holdSubtitle2;dname])
end 
set(ax2,'TickDir','out');
xlim([0 4]);
xticks([]);
ax2.XColor = 'w';
ylim([0 25]);

ax3 = nexttile(3);
hold on
dname = sprintf([feedback ',' ' %0.2f ' char(177) ' %0.2f'], [mean(allDRs(:,25)); std(allDRs(:,25))]);
s3 = swarmchart(ax3, rampNum*ones(1,length(allDRs(:,25))),allDRs(:,25), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c);
title('Recruitment Threshold % EMG MVIC'); holdSubtitle3 = string(ax3.Subtitle.String);
if holdSubtitle3 == ""
    set(ax3.Subtitle,'String', dname)
else
    set(ax3.Subtitle,'String', [holdSubtitle3;dname])
end 
set(ax2,'TickDir','out');
xlim([0 4]); xticks([]); ax2.XColor = 'w'; ylim([0 25]);

ax4 = nexttile(4);
hold on
dname = sprintf([feedback ',' ' %0.2f ' char(177) ' %0.2f'], [mean(allDRs(:,26)); std(allDRs(:,26))]);
s4 = swarmchart(ax4, rampNum*ones(1,length(allDRs(:,26))),allDRs(:,26), 'filled', 'MarkerFaceAlpha', 0.5,'SeriesIndex', c);
title('De-Recruitment Threshold % EMG MVIC'); holdSubtitle4 = string(ax4.Subtitle.String);
if holdSubtitle4 == ""
    set(ax4.Subtitle,'String', dname)
else
    set(ax4.Subtitle,'String', [holdSubtitle4;dname])
end 
set(ax4,'TickDir','out');
xlim([0 4]); xticks([]); ax4.XColor = 'w'; ylim([0 25]);

if keeper1 == 1
    tSave = get(f1Fig, 'Children');
    tSave.Title.String = 'Plantarflexors 20% MVC Ramp 1';
    set(f1Fig, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
    set(f1Fig, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
    saveas(f1Fig, fullfile('G:\Coco Subject Data','Ramp1_Plantarflexors_Deltaf_Summary.pdf'))
end
% end

function [rawdfs, avgdfs, data_MU] = rampSwarm1020(fullFilename)

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
    if size(placeholder,1) < 3
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

%% run random effects model
