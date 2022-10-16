% function Run_BATCH_DeltaF_MT
clearvars % better than clear all

%  [num,txt,raw] = xlsread('Master Subject Data Sheet.xls') ;

% Subject = 'Coco04';
% filenames = cellfun(@num2str,raw,'un',0);

% filenames = dir(fullfile(rootDir, '**', '*TA_*.xls'));

% contraction2Remove = contains({filenames.folder}, 'keep') ;
% filenames(contraction2Remove) = [];

% joint2Remove = ~contains({filenames.folder}, 'Torque Ramps') ;
% filenames(joint2Remove) = [];

%%%% Pick files directory - use wildcard t select specific muscles
% rootDir = cd;
% masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls'); 
% % subject = 'Coco01';
% %subject = 'Coco03_23661_MG_5_v23decomposed';
% [filenames] = getCocoFiles(masterSheet, 'muscle', 'TA', 'feedback', 'Hold', 'contraction', 'DF')

%  joint2Remove = ~contains(filenames, 'Torque Ramps') ;
%  filenames(joint2Remove) = [];

unitMap = colormap;
cutFlicks = 12;
fallingEdge = 12;

[filename,PathName] = uigetfile({'*.mat'},'Select the MATLAB Decomposed file');
fileVariables = {'fsamp','MUPulses', 'JR3Mz','TraceFeedback', 'EMGFeedback'};
load(filename, fileVariables{:});

psf_pdf = fullfile(PathName,'psf_pdf\');

if ~exist(psf_pdf)
    mkdir(psf_pdf);
end

Under = strfind(filename,'_');                      % search for file delimiters

ClusFile = [filename(1:(Under(4)-1)),'.mat'];   % Find cluster file if neccessary
structTAEMG = load(ClusFile, 'EMG');
TAEMG27 = abs(structTAEMG.EMG(27,:)); 
filtTAEMG = movmean(TAEMG27, 500);
EMGFeedback = filtTAEMG-mean(filtTAEMG(1:100));
TAEMG = TAEMG27-mean(TAEMG27(1:100));

solFile = replace(filename, 'TA', 'Sol');
solClusFile = [solFile(1:(Under(4))),'.mat'];   % Find cluster file if neccessary
structSolEMG = load(solClusFile, 'EMGall');
SolEMG27 = abs(structSolEMG.EMGall(27,:));
filtSolEMG = movmean(SolEMG27, 500);
filtSolEMG = filtSolEMG-mean(filtSolEMG(1:100));
SolEMG = SolEMG27-mean(SolEMG27(1:100));

MGFile = replace(filename, 'TA', 'MG');
MGClusFile = [MGFile(1:(Under(4)-1)),'.mat'];   % Find cluster file if neccessary
structMGEMG = load(MGClusFile, 'EMGall');
MGEMG27 = abs(structMGEMG.EMGall(27,:));
filtMGEMG = movmean(MGEMG27, 500); 
filtMGEMG = filtMGEMG-mean(filtMGEMG(1:100));
MGEMG = MGEMG27-mean(MGEMG27(1:100));

MUFiring = SortUnits(MUPulses);
unitSTs = binarySpikeTrain(MUFiring, JR3Mz);
[smoothIDRs,smoothCST,smoothTime] = FilteredUnits(MUFiring, TraceFeedback);

torqueTime = [0:length(JR3Mz)-1]/fsamp;
baseTorque = JR3Mz-mean(JR3Mz(1:fsamp));
if mean(baseTorque) < 0
    realTorque = (baseTorque.*-1)./100;
else
    realTorque = baseTorque./100;
end 

findTorque = realTorque;%(fsamp*10:end-fsamp*10);
findTorque = findTorque<fallingEdge;
[~,torquePeaks] = findpeaks(findTorque*2, 'MinPeakDistance', 10*fsamp); 
Flicks = torqueTime(torquePeaks);
Flicks(Flicks<15 & Flicks>length(realTorque)/fsamp-15);%+10;

periStimTime = [-0.5:0.005:1.5];
stimulusWindow = [-0.5 1.5];

[periTorqueTime,periTorque] = perStimulusAnalog(realTorque, Flicks, torqueTime);

[periTAEMGTime,periTAEMG] = perStimulusAnalog(TAEMG, Flicks, torqueTime);
[periSolEMGTime,periSolEMG] = perStimulusAnalog(SolEMG, Flicks, torqueTime);
[periMGEMGTime,periMGEMG] = perStimulusAnalog(MGEMG, Flicks, torqueTime);

[periEMGFeedbackTime,periEMGFeedback] = perStimulusAnalog(EMGFeedback, Flicks, torqueTime);

[CST_Time, CST_DR, perStimulusCST_Time, perStimulusCST_IDR] = CST_PSF_PSTH(MUFiring, Flicks, stimulusWindow);

realCST_PSTH = histc(perStimulusCST_Time,periStimTime);
figureCST_PSTH = histc(perStimulusCST_Time,periStimTime);

midflick = periTorqueTime>0 == 1 & periTorqueTime <1 == 1;
rFlick = [zeros(size(periTorque,1), 1)];
for i = 1:size(periTorque,1)
    rFlick(i) = mean(periTorque(i,midflick),2);
end 
periTorque(rFlick>cutFlicks,:) = [];

for j = 1:length(MUFiring)
    
    [periStimUnitTime,periStimUnitIDR] = Unit_PSF_PSTH(MUFiring{j},Flicks, stimulusWindow);
    figureUnit_PSTH = histc(periStimUnitTime,periStimTime);
    realUnit_PSTH = histc(periStimUnitTime,periStimTime);

    if j < 20
        unitColor = unitMap(j,:);
        if j == 15 || j == 16
            unitColor = unitMap(j+2,:);
        end
    elseif j > 20
        unitColor = unitMap(j-20,:);
    end
    
    unitFig = figure(j);
    set(unitFig,'Visible','on');
    t = tiledlayout(unitFig, 4,1);% 'TileSpacing','Compact','Padding','Compact');
    title(t,[filename(1:Under(2)), ' Unit ', num2str(j)]);
    
    % Torque Tile
    FlickTQax = nexttile(1);
    hold on
    line([-1 10], [0 0])
    plot(FlickTQax, periTorqueTime, periTorque','.', 'Color',  'k')
    line( [0 0], [-5 10]);
    FlickTQax.XTick = []; FlickTQax.XAxis.Color = 'none'; xlim([-0.5 1.5]); ylim([-10 15]);
    title(FlickTQax, 'Overlayed Torque Triggers');
    
    % EMG Flicks Tile
    EMGsAx = nexttile(2);  
    EMGsAx.YTick = [];
    EMGsAx.XTick = []; EMGsAx.XAxis.Color = 'none';
    plot(EMGsAx, periMGEMGTime, mean(periMGEMG),'Color',  unitMap(5,:), 'LineWidth', 2)
    hold on   
    plot(EMGsAx, periSolEMGTime, mean(abs(periSolEMG)),'Color',  unitMap(7,:), 'LineWidth', 2)
    plot(EMGsAx, periTAEMGTime, mean(abs(periTAEMG)).*.1,'Color', unitMap(1,:), 'LineWidth', 2)
    line([0 0], [-20 50], 'LineWidth', 3, 'Color', 'k'); 
    xlim([-0.5 1.5]); ylim([-20 50]);
    title(EMGsAx, 'Overlayed TA and Sol EMG Triggers');
    
    % PSTH Tile
    PSTHax = nexttile(3)
    bar(PSTHax, periStimTime,(figureCST_PSTH./(sum(realCST_PSTH)/2)).*100, 'FaceColor', [.7 .7 .7], 'EdgeColor', [.7 .7 .7])
    hold on
    bar(PSTHax, periStimTime,(figureUnit_PSTH./sum(realUnit_PSTH)).*100,'FaceColor', unitColor,'EdgeColor', unitColor)
    line( [0 0], [-2 30], 'LineWidth', 3, 'Color', 'k');
    PSTHax.XTick = []; PSTHax.XAxis.Color = 'none';
    ytickformat(PSTHax, 'percentage');
    PSTHax.XLim = [-0.5 1.5]; PSTHax.YLim = [0 1]; yyaxis right ;  
    plot(PSTHax,periStimTime,cumsum(figureCST_PSTH-mean(figureCST_PSTH(periStimTime<0)))./sum(realCST_PSTH),'LineStyle','-','LineWidth', 2, 'Color',  'k')
    plot(PSTHax,periStimTime,cumsum(figureUnit_PSTH-mean(figureUnit_PSTH(periStimTime<-0)))./sum(realUnit_PSTH),'LineStyle','-','LineWidth', 3, 'Color',  unitMap(3,:))
    PSTHax.YAxis(2).Color = 'k';
    
    % PSF Tile
    PSFax = nexttile(4);
    title(PSFax, ' Discharge Rate and CST');
    scatter(PSFax, perStimulusCST_Time, perStimulusCST_IDR, 'filled', 'MarkerFaceColor', [.7 .7 .7])% '.', 'Color',  [.7 .7 .7], 'LineWidth', 10)
    hold on
    scatter(PSFax,periStimUnitTime, periStimUnitIDR, 'filled','MarkerFaceColor', unitColor, 'MarkerEdgeColor', 'k')
    line( [0 0], [-2 30], 'LineWidth', 3, 'Color', 'k');
    ylim([-1 30]); xlim([-0.5 1.5]); ylabel('Discharge Rate'); 
    yyaxis right ; PSFax.YAxis(2).Color = 'k';
     
    [~, timeIndexesUnit] = sort(periStimUnitTime);
    [~, timeIndexesCST] = sort(perStimulusCST_Time);
    
    CorrectedUnitTimes = periStimUnitTime(timeIndexesUnit);
    CorrectedCSTTimes = perStimulusCST_Time(timeIndexesCST);
    
    correctedUnit = periStimUnitIDR(timeIndexesUnit);
    correctedCST = perStimulusCST_IDR(timeIndexesCST);
    
    periUnitSum = cumsum(correctedUnit-mean(correctedUnit(CorrectedUnitTimes<-0)))./sum(correctedUnit);
    periCSTSum = cumsum(correctedCST-mean(correctedCST(CorrectedCSTTimes<-0)))./sum(correctedCST);
    
    plot(PSFax,CorrectedCSTTimes,periCSTSum,'LineWidth', 2, 'Color',  'k')
    plot(PSFax,CorrectedUnitTimes,periUnitSum,'LineStyle','-','LineWidth', 3, 'Color',  unitMap(3,:))

    set(unitFig, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
    set(unitFig, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
    PSFFigname_save = fullfile(PathName,'psf_pdf',[filename(1:end-30),'_Unit_',  num2str(j),'_PSF_PSTH_v2']);
    print (unitFig,'-dpdf',[PSFFigname_save '.pdf']);
    close(unitFig);
end 

    summaryFig = figure(99);
%     set(unitFig,'Visible','off');
    t = tiledlayout(6,2, 'TileSpacing','Compact','Padding','Compact');
    title(t,[filename(1:Under(2)), ' Summary Figure']);
    
    % Torque Tile
    AnalogsAx = nexttile(1);
    AnalogsAx.YTick = [];
    AnalogsAx.XTick = [];
    hold on
    plot(AnalogsAx,torqueTime, realTorque, 'k')
    %plot(AnalogsAx,torqueTime,EMGFeedback,'Color', unitMap(2,:))
    AnalogsAx.XAxis.Color = 'none';
    AnalogsAx.XLim = [0 torqueTime(end)];
    title(AnalogsAx, 'Torque, ');
    
    %Raw EMG Tile
    axEMG2 = nexttile(3,[2,1]);
    hold on
     plot(axEMG2,torqueTime, filtTAEMG.*0.2, 'Color', unitMap(1,:))
    plot(axEMG2,torqueTime, filtMGEMG, 'Color', unitMap(5,:))
    plot(axEMG2,torqueTime, filtSolEMG,  'Color', unitMap(7,:))
   
    legend('TS', 'MG', 'Sol', 'location', 'NorthWest')
    axEMG2.XTick = [];
    axEMG2.XAxis.Color = 'none';
    axEMG2.XLim = [0 length(TAEMG)./fsamp];
    axEMG2.YTick = [];
    axEMG2.YAxis.Color = 'none';
    
    % Flick Tile
    FlickTQax = nexttile(7,[3,1]);
    hold on
    line([-1 2], [0 0])
    plot(FlickTQax, periTorqueTime, periTorque,'.', 'Color',  'k')
    line( [0 0], [-10 5]);
    FlickTQax.XTick = []; FlickTQax.XAxis.Color = 'none'; xlim([-0.5 1.5]); ylim([-5 15]);
    title(FlickTQax, 'Overlayed Torque Triggers');
    
     % Flick Tile
    FlickTQax = nexttile(2);
    hold on
    line([-1 2], [0 0])
    plot(FlickTQax, periTorqueTime, periTorque,'.', 'Color',  'k')
    line( [0 0], [-10 15]);
    FlickTQax.XTick = []; FlickTQax.XAxis.Color = 'none'; xlim([-0.5 1.5]); ylim([-5 15]);
    title(FlickTQax, 'Overlayed Torque Triggers');
    
    % EMG Flicks Tile
    EMGsAx = nexttile(4);
    EMGsAx.YTick = []; EMGsAx.XTick = [];
    plot(EMGsAx, periSolEMGTime, mean(periSolEMG),'Color',  unitMap(7,:))
    hold on
    plot(EMGsAx, periMGEMGTime, mean(periMGEMG),'Color',  unitMap(5,:))
    plot(EMGsAx, periTAEMGTime, mean(periTAEMG).*.1, 'Color', unitMap(1,:))
    line( [0 0], [-10 20]);
    EMGsAx.XTick = []; EMGsAx.XAxis.Color = 'none'; xlim([-0.5 1.5]); ylim([-10 50]);
    title(EMGsAx, 'Overlayed TA and Sol EMG Triggers');
    
    % PSTH Tile
    PSTHax = nexttile(6,[2 1]);
    hold on
    bar(PSTHax, periStimTime,(figureCST_PSTH./(sum(realCST_PSTH))).*100, 'FaceColor', [.7 .7 .7], 'EdgeColor', [.7 .7 .7])
    ytickformat(PSTHax, 'percentage');
    line( [0 0], [-2 30]);
    PSTHax.XLim = [-0.5 1.5]; PSTHax.XTick = []; PSTHax.YLim = [0 1]; yyaxis right ;  
    plot(PSTHax,periStimTime,cumsum(figureCST_PSTH-mean(figureCST_PSTH(periStimTime<-0)))./sum(realCST_PSTH),'LineStyle','-','LineWidth', 2, 'Color',  'k')
    PSTHax.YTick = [];
   
    % PSF Tile
    PSFax = nexttile(10, [2 1]);
    title(PSFax, ' Discharge Rate and CST');
    hold on
    scatter(PSFax, perStimulusCST_Time, perStimulusCST_IDR, 'filled', 'MarkerFaceColor', [.7 .7 .7])% '.', 'Color',  [.7 .7 .7], 'LineWidth', 10)
    line( [0 0], [-2 30]);
    ylim([-1 30]); xlim([-0.5 1.5]); ylabel('Discharge Rate'); yyaxis right ;
    PSFax.YTick= [];
     
    [~, timeIndexesCST] = sort(perStimulusCST_Time);
    CorrectedCSTTimes = perStimulusCST_Time(timeIndexesCST);
    correctedCST = perStimulusCST_IDR(timeIndexesCST);
    periCSTSum = cumsum(correctedCST-mean(correctedCST(CorrectedCSTTimes<-0)))./sum(correctedCST);
    
    plot(PSFax,CorrectedCSTTimes,periCSTSum,'LineWidth', 2, 'Color',  'k')
    
    set(gcf, 'PaperPosition', [0 0 20 15]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [20 15]); %Set the paper to have width 5 and height 5.
    summaryFigname_save = fullfile([filename(1:end-4),'_Summary_PSF_PSTH']);
    print (gcf,'-dpdf',[summaryFigname_save '.pdf']);
