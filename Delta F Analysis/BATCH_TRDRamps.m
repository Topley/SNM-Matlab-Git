clear all
close all
clc

%% Delta F Batch run
clear filenames
% filenames(1).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_20001_kk_6.mat';
% filenames(2).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_20002_kk_6.mat';
 filenames(1).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_70801_kk_5.mat';
% filenames(4).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_24301_kk_5.mat';
% filenames(1).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_40001_kk_5.mat';
% filenames(6).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_40002_kk_5.mat';
direction = 'Forward';
% filenames(2).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_20003_kk_6.mat';
 filenames(2).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_71501_kk_5.mat';
% filenames(9).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_71502_kk_5.mat';
% filenames(2).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_40005_kk_5.mat';
% filenames(11).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_40006_kk_5.mat';
%%
clear filenames
filenames(1).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_20007_kk_6.mat';
filenames(2).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_20008_kk_6.mat';
filenames(3).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_71500_kk_5.mat';
filenames(4).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_24303_kk_5.mat';
filenames(5).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_40007_kk_5.mat';
filenames(6).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_44608_kk_5.mat';

direction = 'Backward';
filenames(7).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_53200_kk_6.mat';
filenames(8).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest05\kk files\TRDTest05_53201_kk_6.mat';
filenames(9).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_73100_kk_5.mat';
filenames(10).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest04\kk files\TRDTest04_73101_kk_5.mat';
filenames(11).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_44609_kk_5.mat';
filenames(12).name = 'F:\Toolkit\Mirror\TRD Hub\TRDTest02b\kk files\TRDTest02b_44610_kk_5.mat';
%% Process Trials

trialRamps = 0;
TimeWindow = [];
Torque = [];
TA_EMG = [];
MG_EMG = [];
Sol_EMG = [];
Auto = [0,0];

for i = 1:size(filenames,2)
    [~, file] = fileparts(filenames(i).name)
    
    [Ramps, AnalogTime, cutAnalog, cutTA, cutMG, cutSol] = TRD_RampsNOTUNITS(filenames(i).name, Auto, direction);
    TimeWindow = [TimeWindow, AnalogTime];
    Torque = [Torque, cutAnalog];
    TA_EMG = [TA_EMG, cutTA];
    MG_EMG = [MG_EMG, cutMG];
    Sol_EMG = [Sol_EMG, cutSol];
    if i == 2
        cols = trialRamps;
    end
    trialRamps = trialRamps + Ramps;
    
    
end

%% Plot Ramps
figure
t = tiledlayout(4,2, 'TileSpacing', 'compact');
nexttile(1)
plot(Torque(:,[1:cols]), 'LineWidth', 1)
hold on
plot(mean(Torque(:,[1:cols]),2), 'k', 'LineWidth', 2.5)
title('Torque',  'FontSize', 15)

nexttile(2)
plot(Torque(:,[cols+1:end]),'LineWidth', 1)
hold on
plot(mean(Torque(:,[cols+1:end]),2), 'k', 'LineWidth', 2.5)
title('Torque',  'FontSize', 15)
hold off


nexttile(3)
plot(Sol_EMG(:,[1:cols]), 'LineWidth', 1)
hold on
plot(mean(Sol_EMG(:,[1:cols]),2), 'k', 'LineWidth', 2.5)
title('Sol',  'FontSize', 15)

nexttile(4)
plot(Sol_EMG(:,[cols+1:end]), 'LineWidth', 1)
hold on
plot(mean(Sol_EMG(:,[cols+1:end]),2), 'k', 'LineWidth', 2.5)
title('Sol',  'FontSize', 15)
hold off

nexttile(5)
plot(MG_EMG(:,[1:cols]), 'LineWidth', 1)
hold on
plot(mean(MG_EMG(:,[1:cols]),2), 'k', 'LineWidth', 2.5)
title('MG',  'FontSize', 15)

nexttile(6)
plot(MG_EMG(:,[cols+1:end]), 'LineWidth', 1)
hold on
plot(mean(MG_EMG(:,[cols+1:end]),2), 'k', 'LineWidth', 2.5)
title('MG',  'FontSize', 15)
hold off

nexttile
plot(TA_EMG(:,[1:cols]), 'LineWidth', 1)
hold on
plot(mean(TA_EMG(:,[1:cols]),2), 'k', 'LineWidth', 2.5)
% ax = gca;
% ax.XTick(1) = [];
% length(ax.XTickLabels);
% time = round(linspace(5, 25, 5));
% ax.XTickLabels = time;
title('TA',  'FontSize', 15)

nexttile
plot(TA_EMG(:,[cols+1:end]), 'LineWidth', 1)
hold on
plot(mean(TA_EMG(:,[cols+1:end]),2), 'k', 'LineWidth', 2.5)
title('TA',  'FontSize', 15)
% ax = gca;
% ax.XTick(1) = [];
% length(ax.XTickLabels);
% time = round(linspace(5, 25, 5));
% ax.XTickLabels = time;
hold off


figAxes = findall(gcf,'type','axes');
for i = 1:8
    figAxes(i).TickLength = [0 0];
    xlim(figAxes(i), [0 25*2048]);
    box(figAxes(i), 'off');
    figAxes(i).FontSize = 15;
    
    if i <7
        ylim(figAxes(i), [0 100]);
    else
        ylim(figAxes(i), [-100 125]);
    end
    
    if i >= 3
        figAxes(i).XTickLabels = [];
    end
    
    if i == 1 || i == 2
        if i == 1
            figAxes(i).YTickLabels = [];
        end
        figAxes(i).XTick(1) = [];
        length(figAxes(i).XTickLabels);
        time = round(linspace(5, 25, 5));
        figAxes(i).XTick = time.* fsamp;
        figAxes(i).XTickLabels = time;
    elseif i == 3 || i == 5 || i == 7
        figAxes(i).YTickLabels = [];
    end
    
end

title(t, [string(direction), " Ramps"])
set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5
%%
saveas(gcf, ['G:\TRD Hub\Test Subject Data\', direction,'_Ramps.png'])
save([direction,'RampsVariables.mat'], 'TimeWindow', 'Torque', 'TA_EMG', 'MG_EMG', 'Sol_EMG', 'cols');

%%
figure(1)
t = tiledlayout(4,2)
figure(2)
t2 = tiledlayout(4,2)
figure(3)
t3 = tiledlayout(1,2)
c = 0;
tSqrE = [];
eSqrE = [];
for jj = 1:4
    
    switch jj
        case 1
            signal = Torque;
        case 2
            signal = Sol_EMG;
        case 3
            signal = MG_EMG;
        case 4
            signal = TA_EMG;
    end
    
    figure(1)
    loops = size(signal,2);
    fitRaw1 = [];
    fitRaw2 = [];
    VarRes1 = [];
    VarRes2 = [];
    rampWindow = [0:length(TimeWindow(:,1))-1] ./ 2048;
    secondbase = length(TimeWindow)./2;
    rawErrors = [];
    nexttile
    for i = 1:cols
        mVar = signal(:,i);
        upStart = secondbase - 8*2048;
        upEnd = secondbase - 2048;
        downEnd = secondbase + 8*2048;
        downStart = secondbase + 2048;
        xUp = rampWindow(upStart:secondbase);
        xDown = rampWindow(secondbase:downEnd);
        m1 = 6;
        x1 = upStart./fsamp;
        y1 = 12;
        targetUp = m1*(xUp - x1) + y1;
        m2 = -6;
        x2 = downStart./fsamp;
        y2 = 54;
        targetDown = m2*(xDown - x2) + y2;
        targetTrace = [targetUp,targetDown(1:end-1)];
        if length(t.Children) > 6
            targetTrace = zeros(length(targetTrace),1)';
        end
        Errors1 = mVar(upStart:downEnd)' - targetTrace;
        rawErrors = [rawErrors,Errors1'];
        RMSE = mean(sqrt(rawErrors(1,:).^2));
        correctTrace = targetTrace - RMSE;
        detrendRes = mVar(upStart:downEnd)' - correctTrace;
        RMSE = mean(sqrt(detrendRes(1,:).^2));
        fitRaw1 = [fitRaw1; RMSE];
        VarRes1 = [VarRes1; detrendRes];
        
        hold on
        plot(rampWindow, mVar, 'Color', [.7 .7 .7])
        plot(rampWindow(upStart:downEnd),targetTrace, 'k', 'LineWidth', 3)
        
    end
    
    nexttile
    for i = cols:loops
        mVar = signal(:,i);
        upStart = secondbase - 8*2048;
        upEnd = secondbase - 2048;
        downEnd = secondbase + 8*2048;
        downStart = secondbase + 2048;
        xUp = rampWindow(upStart:secondbase);
        xDown = rampWindow(secondbase:downEnd);
        m1 = 6;
        x1 = upStart./fsamp;
        y1 = 12;
        targetUp = m1*(xUp - x1) + y1;
        m2 = -6;
        x2 = downStart./fsamp;
        y2 = 54;
        targetDown = m2*(xDown - x2) + y2;
        targetTrace = [targetUp,targetDown(1:end-1)];
        if length(t.Children) > 6
            targetTrace = 5.*ones(length(targetTrace),1)';
        end
        Errors1 = mVar(upStart:downEnd)' - targetTrace;
        rawErrors = [rawErrors,Errors1'];
        RMSE = mean(sqrt(rawErrors(1,:).^2));
        correctTrace = abs(targetTrace - RMSE);
        detrendRes = mVar(upStart:downEnd)' - correctTrace;
        RMSE = mean(sqrt(detrendRes(1,:).^2));
        fitRaw2 = [fitRaw2; RMSE];
        VarRes2 = [VarRes2; detrendRes];
        
        hold on
        plot(rampWindow, mVar, 'Color', [.7 .7 .7])
        plot(rampWindow(upStart:downEnd),targetTrace, 'k', 'LineWidth', 3)
    end
    
    figure(2)
    nexttile
    plot(rampWindow, mean(signal(:,[1:cols]),2), 'k')
    hold on
    sqrEr = mean(sqrt(VarRes1.^2),1);
    plot(rampWindow(upStart:downEnd),detrend(sqrEr), 'SeriesIndex', 3)
    tSqrE = [tSqrE;sqrEr];

    nexttile
    plot(rampWindow, mean(signal(:,[cols:end]),2), 'k')
    hold on
    sqrEr2 = mean(sqrt(VarRes2.^2),1);
    plot(rampWindow(upStart:downEnd),detrend(sqrEr2), 'SeriesIndex', 3)
    eSqrE = [eSqrE;sqrEr2];

    figure(3)
    nexttile(1)
    c = jj * 2 - 2;
    s = swarmchart(jj*2.*ones(length(fitRaw1),1), fitRaw1, 'filled', 'SeriesIndex', c)
    hold on
    s2 = swarmchart(jj*2, mean(sqrEr),'k', 'diamond', 'filled')
    s2.SizeData = 150;
    s.XJitter = 'randn';
    s.XJitterWidth = 0.25;
    xlim([0 10])
    ylim([0 50])
    
    nexttile(2)
    s3 = swarmchart(jj*2.*ones(length(fitRaw2),1), fitRaw2, 'filled', 'SeriesIndex', c)
    hold on
    s4 = swarmchart(jj*2, mean(sqrEr2),'k', 'diamond', 'filled')
    s4.SizeData = 150;
    s3.XJitter = 'randn';
    s3.XJitterWidth = 0.25;
    xlim([0 10])
    ylim([0 50])
    %l = legend(['RMSE ', num2str(mean(sqrEr))]);
    
end
test = mean(tSqrE,2);

legend([[]; num2str(test(1)); []; num2str(test(2)); []; num2str(test(3)); []; num2str(test(4))])
