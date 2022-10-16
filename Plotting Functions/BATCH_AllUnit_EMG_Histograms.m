
% function BATCH_plotMUFiring_AllUnits
filenames = dir('*TA_5.mat')
% filenames = dir('*Sol_5_v23decomposed_MUCLEANED.mat')

for i = 15%1:size(filenames,1)
    try
        plot_AllMUFiring(filenames(i).name)
    catch
        disp(['error with ', filenames(i).name])
    end
    close all
end
% end


function [] = plot_AllMUFiring(toClusterTA)
% plotMUFiring    loads motor units decomposed from cluster; plots units

fsamp = 2048;
[tableau] = Tableau;
under = strfind(toClusterTA, '_');


% clusterFileVariables = {'fsamp','CrossCorrEMG', 'EMG', 'AnalogNames', 'AnalogData'};
toClusterSol = replace(toClusterTA,'TA','Sol');
toClusterMG = replace(toClusterTA,'TA','MG');
% toClusterLG = replace(filenameTA,'TA','LG');

varsTA=  matfile(toClusterTA);
%varsTA=  load(toClusterTA, clusterFileVariables{:});
[~,TA_EMG] = EMGpro(varsTA.EMGall, 'rmsWin',40);


timeEMG = [1:length(TA_EMG)]./fsamp;

varsSol =  matfile(toClusterSol);
%varsSol =  load(toClusterSol, clusterFileVariables{:});
[~,Sol_EMG] = EMGpro(varsSol.EMG, 'rmsWin', 40);

varsMG =  matfile(toClusterMG);
%varsMG =  load(toClusterMG, clusterFileVariables{:});
[~,MG_EMG] = EMGpro(varsMG.EMG, 'rmsWin', 40);

decompTA = [toClusterTA(1:end-4),'_v23decomposed.mat'];
decompSol = replace(decompTA,'TA','Sol');
decompMG = replace(decompTA,'TA','MG');
% toClusterLG = replace(filenameTA,'TA','LG');

otbEMG = varsTA.CrossCorrEMG;
%otbEMG = varsTA.xCorChannel;

analogNamesTA = varsTA.AnalogNames;
analogDataTA = varsTA.AnalogData;

c3dEMGChan = contains(analogNamesTA(:), 'EMG');
c3dEMG = double(analogDataTA(c3dEMGChan,:));

rMxChanC3d = contains(analogNamesTA(:), '.rMx');
rMx = double(analogDataTA(rMxChanC3d,:));
rMx = rMx-round(mean(rMx(1:500)));

rMyChanC3d = contains(analogNamesTA(:), '.rMy');
rMy = double(analogDataTA(rMyChanC3d,:));
rMy = rMy-round(mean(rMy(1:500)));

rFxChanC3d = contains(analogNamesTA(:), '.rFx');
rFx = double(analogDataTA(rFxChanC3d,:));
rFx = rFx-round(mean(rFx(1:500)));

rFyChanC3d = contains(analogNamesTA(:), '.rFy');
rFy = double(analogDataTA(rFyChanC3d,:));
rFy = rFy-round(mean(rFy(1:500)));

rFzChanC3d = contains(analogNamesTA(:), '.rFz');
rFz = double(analogDataTA(rFzChanC3d,:));
rFz = rFz-round(mean(rFz(1:500)));

rCOPx = (rMy + rFx)./(rFz + 4);
rCOPy = (rMx + rFy)./(rFz + 4);
c3dTime = [1:length(rFx)]./960;

[fsampC3d,difT] = emg_sync(otbEMG, c3dEMG, fsamp, 960, 0);
if difT < 0
    difT = difT*-1;
end

try
    %TAPulses = load(decompTA, 'MUPulses');
    TAPulses = matfile(decompTA);
    TAFiring = SortUnits(TAPulses.MUPulses);
    for i = 1:length(TAFiring)
        delIndx = (TAFiring{i} <= difT*fsampC3d);
        if max(delIndx) == 1
        TAFiring{i}(delIndx) = [];
        TAFiring{i} = TAFiring{i}-TAFiring{i}(1);
        end 
    end
end

try
    %SolPulses = load(decompSol, 'MUPulses');
    SolPulses = matfile(decompSol);
    SolFiring = SortUnits(SolPulses.MUPulses);
    for j = 1:length(SolFiring)
         delIndx = (SolFiring{j} <= difT*fsampC3d);
        if max(delIndx) == 1
        SolFiring{j}(delIndx) = [];
        SolFiring{j} = SolFiring{j}-SolFiring{j}(1);
        end 
    end
end

try
    %MGPulses = load(decompMG, 'MUPulses');
    MGPulses = matfile(decompMG);
    MGFiring = SortUnits(MGPulses.MUPulses);
     for ii = 1:length(MGFiring)
         delIndx = (MGFiring{ii} <= difT*fsampC3d);
        if max(delIndx) == 1
        MGFiring{ii}(delIndx) = [];
        MGFiring{ii} = MGFiring{ii}-MGFiring{ii}(1);
        end
    end
end

offsetT = find(timeEMG <= difT);
timeEMG(offsetT) = [];
timeEMG = timeEMG-timeEMG(1);

TA_EMG(offsetT) = [];
Sol_EMG(offsetT) = [];
MG_EMG(offsetT) = [];

h=figure('visible', 'off');
t = tiledlayout(15, 1, 'TileSpacing', 'tight');
t1 = nexttile([5,1]);
t1.YAxis.TickLabels = []; t1.XAxis.Color = 'none'; t1.YAxis.Color = 'none';
 
% ax2 = axes(t); ax2.Layout.Tile = 2;
% CST = sort(cell2mat(TAFiring));
% scatter(ax2,CST(2:end)./fsamp, 1./diff(CST).*fsamp, 'filled', 'SeriesIndex', 1)

ax3 = axes(t); ax3.Layout.Tile = 1; ax3.Layout.TileSpan = [4 1];
try
    for i = 1:size(TAFiring,2)
        tsp = TAFiring{i}/fsamp;
        spikes = numel(tsp);
        for j = 1:spikes
            hold on
            line(ax3,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(2,:))
        end
    end
    ylim(ax3, [-i-1 0]);
end
ax1 = axes(t); ax1.Layout.Tile = 4; ax1.Layout.TileSpan = [2 1];
plot(ax1,timeEMG, TA_EMG, 'SeriesIndex', 1)
ylim(ax1, [0 max(TA_EMG)]);

t2 = nexttile([5,1]);
t2.YAxis.TickLabels = []; t2.XAxis.Color = 'none'; t2.YAxis.Color = 'none';
% ax5 = axes(t); ax5.Layout.Tile = 7;
% CST = sort(cell2mat(SolFiring));
% scatter(ax5,CST(2:end)./fsamp, diff(CST), 'filled', 'SeriesIndex', 7)

ax6 = axes(t); ax6.Layout.Tile = 6; ax6.Layout.TileSpan = [4 1];
try
    for i = 1:size(SolFiring,2)
        tsp = SolFiring{i}/fsamp;
        spikes = numel(tsp);
        for j = 1:spikes
            hold on
            line(ax6,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(8,:))
        end
    end
    ylim(ax6, [-i-1 0]);
end
ax4 = axes(t); ax4.Layout.Tile =9; ax4.Layout.TileSpan = [2 1];
plot(ax4,timeEMG, Sol_EMG,'SeriesIndex',7)
ylim(ax4, [0 max(Sol_EMG)]);

t3 = nexttile([5,1]);
t3.YAxis.TickLabels = []; t3.YAxis.Color = 'none';
% ax8 = axes(t); ax8.Layout.Tile = 11;
% CSTm = sort(cell2mat(MGFiring));
% scatter(ax8,CST(2:end)./fsamp, diff(CST), 'filled', 'SeriesIndex', 5)

ax9 = axes(t); ax9.Layout.Tile = 11; ax9.Layout.TileSpan = [4 1];
try
    for i = 1:size(MGFiring,2)
        tsp = MGFiring{i}/fsamp;
        spikes = numel(tsp);
        for j = 1:spikes
            line(ax9,[tsp(j) tsp(j)], [-i -i-1], 'Color', tableau(6,:))
        end
    end
ylim(ax9, [-i-1 0]);    
end
ax7 = axes(t); ax7.Layout.Tile = 14; ax7.Layout.TileSpan = [2 1];
plot(ax7,timeEMG, MG_EMG,'SeriesIndex',5)
ylim(ax7, [0 max(MG_EMG)]);
xlabel(t,'Time');

topAx = nexttile(t, 'North');
hold on
plot(c3dTime, rMx + mean(rMx(1:fsamp))*2, 'k')
plot(c3dTime, rMy, 'Color', [.5 .5 .5])

figAxes = findall(gcf,'type','axes');
for i = 1:length(figAxes)
    xlim(figAxes(i), [0 length(TA_EMG)./fsamp])
    set(figAxes(i), 'Color', 'none');
    set(figAxes(i).XAxis, 'Color', 'none');
    set(figAxes(i).YAxis, 'Color', 'none');
    set(figAxes(i), 'Color', 'none');
end
set(t1,'Color', 'w');set(t2,'Color', 'w');set(t3,'Color', 'w');set(t3.XAxis,'Color', 'k'); set(topAx,'Color', 'w');
% set(h, 'Visible', 'on')

title([toClusterTA(1:(under(2))),'AllUnits_v23decomposed.mat'])

set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5
savefilename = strcat([toClusterTA(1:(under(2))),'AllUnits_v23decomposed'], '_CHECK.pdf');
print (h,'-dpdf',savefilename);
close(h)
%
% if check == 1
%     if isfile(strcat(filenameTA(1:end-4), '_MUCLEANED_REVIEW.pdf')) == 1;
%         close(h)
%         return
%     end
%     savefilename = strcat(filenameTA(1:end-4), '_CHECK.pdf');
%     print (h,'-dpdf',savefilename);
%     pause(1)
%     close(h)
%
% elseif check == 0
%     try
%         warning('off', 'MATLAB:DELETE:FileNotFound');
%         delete(strcat(filenameTA(1:end-4),'CHECK.pdf'));
%     catch
%     end
%     %      if isfile(strcat(filename(1:end-4), '_REVIEW.pdf')) == 1;
%     %         close(h)
%     %         return
%     %     end
%     savefilename = strcat(filenameTA(1:end-4), '_REVIEW.pdf');
%     print (h,'-dpdf',savefilename);
%     pause(1)
%     close(h)
% else
% end

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
