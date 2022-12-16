clear all
close all
clc
%% All Muscle Summary Raster Plots

% Must be in outermost directory containing all study subject data
SubjectID = 'TRDTest05';
c = Tableau;
SubjectDir = cd;
rootDir = fullfile(SubjectDir, SubjectID);
kkFiles = dir(fullfile(rootDir,'kk files', '*.mat'));

for ii=1:length(kkFiles)
    % plotMUFiring    loads motor units from all decomposed files and variables
    % from toCluster files; plots raster of units, IDRs, average EMG, and some Treadmill Forces
    
    kkFile = fullfile(kkFiles(ii).folder, kkFiles(ii).name);
    [~, clusterFile, decompFile, cleanFile] = Get_TrialFiles(kkFile);
    
    [~, TrialID] = fileparts(kkFile);
    TrialID = TrialID(1:end-5);
    [fileDir, ~]  = fileparts(decompFile);
    fsamp = 2048;
    kkVars = matfile(kkFile);
    copS = kkVars.COP;
    cpML = copS.WeightedX;
    cpAP = copS.WeightedY;
    fpS = kkVars.ForcePlates;
    RightFz = fpS.RightFz;
    TimeDelay = kkVars.TimeDelay;
    TimeStop = kkVars.TimeStop;
    
    try
        Markers = kkVars.Markers;
        Sternum = -1.*Markers.STR(:,2);
        Sternum = Sternum(1:end-1);
    end
    
    unpadAP = RightFz ~= 0;
    paddAP = RightFz == 0;
    offsetR = min(RightFz(unpadAP));
    RightFz(paddAP) = offsetR;
    rFz = RightFz;
    RightFz = fpS.RightFz;
    
    unpadAP = cpAP ~= 0;
    paddAP = cpAP == 0;
    offset = min(cpAP(unpadAP));
    cpAP(paddAP) = offset;
    copAP = cpAP;
    unpadML = cpML ~= 0;
    paddML = cpML == 0;
    offset2 = min(cpML(unpadML));
    cpML(paddML) = offset2;
    copML = cpML;
    
    clusterTA = replace(clusterFile, 'kk', 'TA');
    varsTA = matfile(clusterTA);
    clusterSol = replace(clusterFile, 'kk', 'Sol');
    varsSol = matfile(clusterSol);
    clusterMG = replace(clusterFile, 'kk', 'MG');
    varsMG = matfile(clusterMG);
    
    try
        decomposedTA = replace(decompFile, 'kk', 'TA');
    end
    
    try
        decomposedSol = replace(decompFile, 'kk', 'Sol');
    end
    
    try
        decomposedMG = replace(decompFile, 'kk', 'MG');
    end
    cleanTA = replace(cleanFile, 'kk', 'TA');
    cleanSol = replace(cleanFile, 'kk', 'Sol');
    cleanMG = replace(cleanFile, 'kk', 'MG');
    
    try
        unitsTA = matfile(cleanTA);
        TAPulses = unitsTA.MUPulses;
        TAFiring = SortUnits(TAPulses);
    catch
        try
            unitsTA = matfile(decomposedTA);
            TAPulses = unitsTA.MUPulses;
            TAFiring = SortUnits(TAPulses);
        end
    end
    
    try
        unitsSol = matfile(cleanSol);
        SolPulses = unitsSol.MUPulses;
        SolFiring = SortUnits(SolPulses);
    catch
        try
            unitsSol = matfile(decomposedSol);
            SolPulses = unitsSol.MUPulses;
            SolFiring = SortUnits(SolPulses);
        end
    end
    
    try
        unitsMG = matfile(cleanMG);
        MGPulses = unitsMG.MUPulses;
        MGFiring = SortUnits(MGPulses);
    catch
        try
            unitsMG = matfile(decomposedMG);
            MGPulses = unitsMG.MUPulses;
            MGFiring = SortUnits(MGPulses);
        end
    end
    
    
    Sol_EMG = varsSol.EMG;
    Sol_EMG = Sol_EMG - mean(Sol_EMG, 2);
    SOL_rectified=abs(Sol_EMG);
    Sol_EMG=mean(SOL_rectified,1);
    
    MG_EMG = varsMG.EMG;
    MG_rectified=abs(MG_EMG);
    MG_EMG=mean(MG_rectified,1);
    MG_EMG = MG_EMG - mean(MG_EMG, 2);
    
    TA_EMG = varsTA.EMG;
    TA_EMG = TA_EMG - mean(TA_EMG, 2);
    TA_rectified=abs(TA_EMG);
    TA_EMG=mean(TA_rectified,1);
    timeEMG = [0:length(copAP)-1]./fsamp;
    copTime = timeEMG;
    
    h=figure('visible', 'on');
    t = tiledlayout(18, 1, 'TileSpacing', 'tight');
    t1 = nexttile([3,1]); t1.YAxis.TickLabels = []; t1.XAxis.Color = 'none'; t1.YAxis.Color = 'none';
    t2 = nexttile([5,1]); t2.YAxis.TickLabels = []; t2.XAxis.Color = 'none'; t2.YAxis.Color = 'none';
    t3 = nexttile([5,1]); t3.YAxis.TickLabels = []; t3.XAxis.Color = 'none'; t3.YAxis.Color = 'none';
    t4 = nexttile([5,1]); t4.YAxis.TickLabels = []; t4.YAxis.Color = 'none';
    
    ax1 = axes(t); ax1.Layout.Tile = 1; ax1.Layout.TileSpan = [3 1];
    z = zeros(size(copAP))';
    surface([timeEMG;timeEMG],[copAP';copAP'],[z;z],[timeEMG;timeEMG],...
        'facecolor','none','edgecolor','interp','linewidth',1);
    colormap(ax1,parula);
    
    hold(ax1, 'on');
    plot(ax1, copTime, Sternum./10, 'k')
    
    ax3 = axes(t); ax3.Layout.Tile = 4; ax3.Layout.TileSpan = [4 1];
    ax4 = axes(t); ax4.Layout.Tile = 7; ax4.Layout.TileSpan = [2 1];
    ax5 = axes(t); ax5.Layout.Tile = 9; ax5.Layout.TileSpan = [4 1];
    ax6 = axes(t); ax6.Layout.Tile = 12; ax6.Layout.TileSpan = [2 1];
    ax7 = axes(t); ax7.Layout.Tile = 14; ax7.Layout.TileSpan = [4 1];
    ax8 = axes(t); ax8.Layout.Tile = 17; ax8.Layout.TileSpan = [2 1];
    
    l1 = plot(ax4,timeEMG, TA_EMG, 'SeriesIndex', 7)
    emgC1 = get(l1, 'Color');
    l1.Color = [emgC1(1) emgC1(2) emgC1(3) .5];
    ylim(ax4, [0 max(TA_EMG)]);
    
    l2 = plot(ax6,timeEMG, MG_EMG, 'SeriesIndex', 9)
    emgC2 = get(l2, 'Color');
    l2.Color = [emgC2(1) emgC2(2) emgC2(3) .5];
    ylim(ax6, [0 max(MG_EMG)]);
    
    l3 = plot(ax8,timeEMG, Sol_EMG, 'SeriesIndex',1)
    emgC3 = get(l3, 'Color');
    l3.Color = [emgC3(1) emgC3(2) emgC3(3) .5];
    ylim(ax8, [0 max(Sol_EMG)]);
    
    try
        for i = 1:size(TAFiring,2)
            tsp = TAFiring{i}/fsamp;
            spikes = numel(tsp);
            tsp(tsp < 0) = [];
            for j = 1:spikes
                hold on
                line(ax3,[tsp(j) tsp(j)], [-i -i-1], 'Color', c(8,:))
            end
        end
        ylim(ax3, [-i-2 -1]);
    end
    
    try
        for i = 1:size(SolFiring,2)
            tsp = SolFiring{i}/fsamp;
            tsp(tsp < 0) = [];
            spikes = numel(tsp);
            for j = 1:spikes
                hold on
                line(ax7,[tsp(j) tsp(j)], [-i -i-1], 'Color',  c(2,:))
            end
        end
        ylim(ax7, [-i-2 -1]);
    end
    
    try
        for i = 1:size(MGFiring,2)
            tsp = MGFiring{i}/fsamp;
            tsp(tsp < 0) = [];
            spikes = numel(tsp);
            for j = 1:spikes
                line(ax5,[tsp(j) tsp(j)], [-i -i-1], 'Color', c(10,:))
            end
        end
        ylim(ax5, [-i-2 -1]);
    end
    xlabel(t,'Time');
    
    figAxes = findall(gcf,'type','axes');
    for i = 1:length(figAxes)
        xlim(figAxes(i), [TimeDelay TimeStop])
        set(figAxes(i), 'Color', 'none');
        set(figAxes(i).XAxis, 'Color', 'none');
        set(figAxes(i).YAxis, 'Color', 'none');
        set(figAxes(i), 'Color', 'none');
    end
    set(t1,'Color', 'w');set(t2,'Color', 'w'); set(t3,'Color', 'w'); set(t4,'Color', 'w');set(t4.XAxis,'Color', 'k');
    
    %set(h, 'Visible', 'on')
    title(t, TrialID);
    set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5
    savefilename = fullfile(fileDir,strcat([TrialID,'_AllUnits_Summary'], '.pdf'));
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
