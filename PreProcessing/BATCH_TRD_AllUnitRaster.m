clear all
close all
rootDir = cd;
filenames = dir(fullfile(rootDir, 'to_cluster', '*TA_5.mat'));
% filenames = dir('*Sol_5_v23decomposed_MUCLEANED.mat')
tic
for i = 1:size(filenames,1)
    
    try
        disp(['Processing ', filenames(i).name])
        plot_AllMUFiring(fullfile(filenames(i).folder,filenames(i).name))
    catch
        disp(['   Cant process: ', filenames(i).name])
    end
    
end
toc

function plot_AllMUFiring(toClusterTA)
% plotMUFiring    loads motor units from all decomposed files and variables
% from toCluster files; plots raster of units, IDRs, average EMG, and some Treadmill Forces

[fileDirectory, trialID, mucleanFile, decompFile, clusterFile, kkFile, ~] = Set_Directories(toClusterTA, 'TRD');
[~,filename] = fileparts(clusterFile);
fileDir  = [fileDirectory,'\decomposed'];
[tableau] = Tableau;

under = find(toClusterTA == '_');

toClusterSol = replace(toClusterTA,'TA','Sol');
toClusterMG = replace(toClusterTA,'TA','MG');

varsTA =  matfile(toClusterTA);
[~,TA_EMG] = EMGpro(varsTA.EMGall, 'channel', 27, 'filter', {'filtRMS', 500, 500});
fsamp = varsTA.fsamp;
timeEMG = [0:length(TA_EMG)-1]./fsamp;

varsSol =  matfile(toClusterSol);
[~,Sol_EMG] = EMGpro(varsSol.EMGall, 'channel', 27, 'filter', {'filtRMS', 500, 500});

varsMG =  matfile(toClusterMG);
[~,MG_EMG] = EMGpro(varsMG.EMGall, 'channel', 27, 'filter', {'filtRMS', 500, 500});

decompTA = decompFile;
decompSol = replace(decompTA,'TA','Sol');
decompMG = replace(decompTA,'TA','MG');

varsKK =  matfile(kkFile);
varsKK = varsKK.c3d_Data;
ForcePlates = varsKK.ForcePlates;
COP = varsKK.COP;
copAP = COP.WeightedY;
copML = COP.WeightedX;
copTime = [0:length(copAP)-1]./fsamp;

% timeDiffTA = (length(TA_EMG) - length(copAP))/fsamp;
% timeDiffSol = (length(Sol_EMG) - length(copAP))/fsamp;
% timeDiffMG = (length(MG_EMG) - length(copAP))/fsamp;
% 
% timeTA = timeEMG - timeDiffTA;
% timeSol = timeEMG - timeDiffSol;
% timeMG = timeEMG - timeDiffMG;
% 
% TA_EMG(timeTA < 0) = [];
% Sol_EMG(timeSol < 0) = [];
% MG_EMG(timeMG < 0) = [];
% 
% shiftEMG = timeEMG - timeDiffTA;
% shiftEMG(shiftEMG< 0) = [];
% timeEMG = shiftEMG; 

try
    TAPulses = matfile(decompTA);
    TAFiring = SortUnits(TAPulses.MUPulses);
end

try
    SolPulses = matfile(decompSol);
    SolFiring = SortUnits(SolPulses.MUPulses);
end

try
    MGPulses = matfile(decompMG);
    MGFiring = SortUnits(MGPulses.MUPulses);
end

h=figure('visible', 'off');
t = tiledlayout(17, 1, 'TileSpacing', 'tight');
t1 = nexttile([2,1]); t1.YAxis.TickLabels = []; t1.XAxis.Color = 'none'; t1.YAxis.Color = 'none';
t2 = nexttile([5,1]); t2.YAxis.TickLabels = []; t2.XAxis.Color = 'none'; t2.YAxis.Color = 'none';
t3 = nexttile([5,1]); t3.YAxis.TickLabels = []; t3.XAxis.Color = 'none'; t3.YAxis.Color = 'none';
t4 = nexttile([5,1]); t4.YAxis.TickLabels = []; t4.YAxis.Color = 'none';

ax1 = axes(t); ax1.Layout.Tile = 1;
ax2 = axes(t); ax2.Layout.Tile = 2;
plot(ax1,copTime, copAP, 'k')
plot(ax2,copTime, copML./100, 'Color', [.5 .5 .5])

ax3 = axes(t); ax3.Layout.Tile = 3; ax3.Layout.TileSpan = [4 1];
ax4 = axes(t); ax4.Layout.Tile = 6; ax4.Layout.TileSpan = [2 1];
ax5 = axes(t); ax5.Layout.Tile = 8; ax5.Layout.TileSpan = [4 1];
ax6 = axes(t); ax6.Layout.Tile = 11; ax6.Layout.TileSpan = [2 1];
ax7 = axes(t); ax7.Layout.Tile = 13; ax7.Layout.TileSpan = [4 1];
ax8 = axes(t); ax8.Layout.Tile = 16; ax8.Layout.TileSpan = [2 1];

plot(ax4,timeEMG, TA_EMG, 'SeriesIndex', 1)
ylim(ax4, [0 max(TA_EMG)]);
plot(ax6,timeEMG, Sol_EMG,'SeriesIndex',7)
ylim(ax6, [0 max(Sol_EMG)]);

plot(ax8,timeEMG, MG_EMG,'SeriesIndex',5)
ylim(ax8, [0 max(MG_EMG)]);

try
    for i = 1:size(TAFiring,2)
        tsp = TAFiring{i}/fsamp;
        spikes = numel(tsp);
        tsp(tsp < 0) = [];
        for j = 1:spikes
            hold on
            line(ax3,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(2,:))
        end
    end
    ylim(ax3, [-i-1 0]);
end

try
    for i = 1:size(SolFiring,2)
        tsp = SolFiring{i}/fsamp;
        tsp(tsp < 0) = [];
        spikes = numel(tsp);
        for j = 1:spikes
            hold on
            line(ax5,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(8,:))
        end
    end
    ylim(ax5, [-i-1 0]);
end

try
    for i = 1:size(MGFiring,2)
        tsp = MGFiring{i}/fsamp;
        tsp(tsp < 0) = [];
        spikes = numel(tsp);
        for j = 1:spikes
            line(ax7,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(6,:))
        end
    end
    ylim(ax7, [-i-1 0]);
end
xlabel(t,'Time');

figAxes = findall(gcf,'type','axes');
for i = 1:length(figAxes)
    xlim(figAxes(i), [0 length(TA_EMG)./fsamp])
    set(figAxes(i), 'Color', 'none');
    set(figAxes(i).XAxis, 'Color', 'none');
    set(figAxes(i).YAxis, 'Color', 'none');
    set(figAxes(i), 'Color', 'none');
end
set(t1,'Color', 'w');set(t2,'Color', 'w'); set(t3,'Color', 'w'); set(t4,'Color', 'w');set(t4.XAxis,'Color', 'k');
 set(h, 'Visible', 'on')
title(t, [trialID,'_AllUnits_v23decomposed.mat'])
set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5
savefilename = fullfile(fileDir,strcat([trialID,'_AllUnits_v23decomposed'], '_CHECK.pdf'));
print (h,'-dpdf',savefilename);
close(h)
end


function [tableau] = Tableau
tableau =     [0.1216    0.4667    0.7059
    0.6824    0.7804    0.9098
    1.0000    0.4980    0.0549
    1.0000    0.7333    0.4706
    0.1725    0.6275    0.1725
    0.5961    0.8745    0.5412
    0.8392    0.1529    0.1569
    1.0000    0.5961    0.5882
    0.5804    0.4039    0.7412
    0.7725    0.6902    0.8353
    0.5490    0.3373    0.2941
    0.7686    0.6118    0.5804
    0.8902    0.4667    0.7608
    0.9686    0.7137    0.8235
    0.4980    0.4980    0.4980
    0.7804    0.7804    0.7804
    0.7373    0.7412    0.1333
    0.8588    0.8588    0.5529
    0.0902    0.7451    0.8118
    0.6196    0.8549    0.8980]; % load colormap
end
