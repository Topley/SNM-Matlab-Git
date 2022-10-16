clearvars('-except', 'f1Fig', 'ax1') %%close all
fsamp = 2048;
LW = 1;
figColor = 5;
feedback = 'Coco';
Muscle1 = 'TA';
Muscle2 = 'MG';

rootDir = cd;
filenames = dir(fullfile(rootDir, '**', '*ALLUnits_pooled_coher*.mat'));

trials2Remove = contains({filenames.folder}, 'Test Subjects') ;
filenames(trials2Remove) = [];

masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls');
% subject = 'Coco01';
%subject = 'Coco03_23661_MG_5_v23decomposed';
[files] = getCocoFiles(masterSheet, 'muscle', 'TA', 'feedback', 'Hold', 'contraction', 'DF')
trimFiles = cellfun(@(x) x(1:12),files,'UniformOutput',false);
trials2Remove = ~contains({filenames.name}, trimFiles) ;
filenames(trials2Remove) = [];

avgZ = [];
avgSolZ = [];
avgMGZ =[];
TAexplained = [];
MGexplained = [];
Varexplained = [];
TASolexplained = [];
TAMGexplained = [];

if exist('f1Fig') == 0
    f1Fig = figure(1);%('Visible', 'off');
    t = tiledlayout(1,1);
end
f1Fig = figure(1)
for i = 1:length(filenames)
    try
        load(fullfile(filenames(i).folder,filenames(i).name), 'data', 'fsamp')
        ax1 = nexttile(1)
        hold(ax1,'on')
        F = data.TAF;
        Z = data.TAZ;
        plot(ax1,F,Z, 'SeriesIndex', 4)
        avgZ = [avgZ,Z];
        Varexplained = [Varexplained,data.TAexplained(1:3,1)];
    end
end

try
    avgZp = mean(avgZ,2);
    plot(ax1,F,avgZp,'SeriesIndex', 3, 'LineWidth', 2)
end

% delta = trapz(F(1:50), avgZ(1:50));    % 0 - 5Hz
% alpha = trapz(F(50:150), avgZ(50:150));    % 5 - 15Hz
% beta = trapz(F(150:350), avgZ(150:350));    % 15 - 35Hz
xlim(ax1,[0 40])
ylim(ax1,[0 1])
title([feedback,' ',  Muscle1, '-', Muscle2]);
% subtitle({['Delta ' num2str(mean(delta)) ' Alpha ' num2str(mean(alpha)) ' Beta ' num2str(mean(beta))];[' Variance Explained fPC ', num2str(mean(Varexplained(1,:)))]})
mean(Varexplained(2,:))
% NoiseCL = max(avgTAZ(TAF>100));
% CL = abs(1- (1-0.95)^(1/(data.TAMULength)/(fsamp*LW)-1))/10241;
% plot(ax1,TAF,ones(1,length(TAF))*CL,'r', 'LineWidth', 2);
% plot(ax1,TAF,ones(1,length(TAF))*NoiseCL,'k', 'LineWidth', 2);
% yticks(ax1,[0 0.5 1]);

% avgSolZ = mean(avgSolZ,2);
% plot(ax5,SolF,avgSolZ,'SeriesIndex', 7, 'LineWidth', 2)
% yticks(ax5,[0 0.5 1]);
%
% avgpoolZ = mean(avgpoolZ,2);
% plot(ax3,TAF,avgpoolZ,'Color', 'k', 'LineWidth', 2)
% yticks(ax3,[0 0.5 1]);
%
% %xlabel(ax,'Frequency (Hz)')
% %ylabel(ax,'Coherence (z-transform)')
% %title(ax,[Muscle,'-', Muscle,' Group Pooled Coherence'])
%
% meanAlpha = mean(alpha);
% meanDelta = mean(delta);
% meanBeta = mean(beta);
% % hold(ax,'off');
% set(f1Fig,'Visible','on');

%%
tSave = get(f1Fig, 'Children');
%tSave.Title.String = 'All Subjects Pooled Coherence';
set(f1Fig, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(f1Fig, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(f1Fig, fullfile('G:\Coco Subject Data\Group Figures\','TAMG_Coco_IntraMuscularCoherence.pdf'))

