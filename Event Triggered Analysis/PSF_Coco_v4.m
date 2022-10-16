PSF_folder = 'F:\Toolkit\Mirror\CoCo09B\Decomposed\PSFs';

[filename,PathName] = uigetfile({'*.mat'},'Select the MATLAB Decomposed file');
load(filename);

Under = strfind(filename,'_');                      % search for file delimiters
solFile = replace(filename, 'TA', 'Sol');
solClusFile = [solFile(1:(Under(4))),'.mat'];   % Find cluster file if neccessary
structSolEMG = load(solClusFile, 'EMG');
[filtSolEMG] = allChanRmsEMG(structSolEMG.EMG);
% filtSolEMG(filtSolEMG>10) = mean(filtSolEMG);%  = mean(filtMGEMG);

MGFile = replace(filename, 'TA', 'MG');
MGClusFile = [MGFile(1:(Under(4)-1)),'.mat'];   % Find cluster file if neccessary
structMGEMG = load(MGClusFile, 'EMG');
[filtMGEMG] = allChanRmsEMG(structMGEMG.EMG);

% filtMGEMG(filtMGEMG>10) = mean(filtMGEMG);%  = mean(filtMGEMG);
% if max(filtAntagonistEMG)>20
% structAntagonistEMG = load(antagonistClusFile, 'EMGall');
% antagonistEMG = structAntagonistEMG.EMGall(27,:);
% [filtAntagonistEMG] = vectorFiltEMG(antagonistEMG);
% end 

[filtEMGall] = allChanRmsEMG(EMG);

MUFiring = SortUnits(MUPulses);

unitSTs = binarySpikeTrain(MUFiring);
% smoothTime = [1:length(unitSTs)]/fsamp;
% smoothDRs = fftfilt(hanning(fsamp),unitSTs')'.*2;

[smoothIDRs,smoothCST,smoothTime] = FilteredUnits(MUFiring, TraceFeedback);

torqTime = [1:length(TraceFeedback)]/fsamp;

traceBase = TraceFeedback-mean(TraceFeedback(1:fsamp));
[B,A] = butter(2,1/(fsamp/2));
trace_butter = filtfilt(B,A,traceBase')';
traceReal = trace_butter./3+6;
traceMax = traceReal.*-1+6;

realEMG = EMGFeedback+abs(min(EMGFeedback));

torqueBase = JR3Mz-mean(JR3Mz(1:3*fsamp));
torque_butter = filtfilt(B,A,torqueBase')';
realTorque = torque_butter./3+4;
timeTQ = [0:length(traceReal)-1]/fsamp;

[periTorqueTime,periTorqueSignal] = perStimulusAnalog(realTorque, Flicks);
[periTAEMGTime,periTAEMG] = perStimulusAnalog(realEMG, Flicks);
[periSolEMGTime,periSolEMG] = perStimulusAnalog(filtSolEMG, Flicks);
[periMGEMGTime,periMGEMG] = perStimulusAnalog(filtMGEMG, Flicks);

TF = islocalmax(traceMax, 'MinSeparation', 5*fsamp, 'MinProminence', 2);
Flicks = timeTQ(TF);
Flicks(Flicks<5) = [];
Flicks(Flicks>210) = [];
Flicks = Flicks-1.25;
realTorque = (TorqueFeedback-5).*2+6.5;

[periTorqueTime,periTorque] = perStimulusAnalog(realTorque, Flicks);
[periTAEMGTime,periTAEMG] = perStimulusAnalog(realEMG, Flicks);
[periSolEMGTime,periSolEMG] = perStimulusAnalog(filtSolEMG, Flicks);
[periMGEMGTime,periMGEMG] = perStimulusAnalog(filtMGEMG, Flicks);
[periCSTTime,periCST] = perStimulusAnalog(smoothCST, Flicks);
[CST_Stim_Time, CST_Stim_DR, CST_PSTH, CST_Timebase] = CST_PSF_PSTH(MUFiring, Flicks);

%h = histfit(CST_Stim_DR, length(CST_Stim_DR), 'Kernel')

for j = 1:length(MUFiring)
    
    [Unit_Stim_Time,Unit_Stim_IDR,Unit_Timebase,Unit_PSTH] = Unit_PSF_PSTH(MUFiring{j},Flicks);
    avgDRs = [];
    unitFig = figure(j);
    unitColor = colormap;
%     set(unitFig,'Visible','off');
    t = tiledlayout(3,2);
    title(t,[filename(1:Under(2)), ' Unit ', num2str(j)]);
    
    % Torque Tile
    Torqax = nexttile(1);
    hold on
    plot(timeTQ, realTorque, 'k')
    xlim([0 length(timeTQ)/fsamp]);
    ylim([-10 20]);
    ylabel('Torque');
    title(Torqax, 'Torque');
    
    % Flick Tile
    Flickax = nexttile(2);
%     loopTorque = [];
%     for i = 1:length(Flicks)
%         loopTQ = [];
%         loopTQ(1,:) = realTorque(timeTQ>Flicks(i)-2 & timeTQ<Flicks(i)+2);
%         loopTQ(2,:) = timeTQ(timeTQ>Flicks(i)-2 & timeTQ<Flicks(i)+2) - Flicks(i);
%         loopTorque = [loopTQ loopTorque];
%     end
%     
%     tq_Stim_Time = loopTorque(2,:);
%     tq_Stim_DR = loopTorque(1,:);
%    plot(Flickax, tq_Stim_Time, tq_Stim_DR,'.', 'Color',  [.5 .5 .5])
    plot(Flickax, periTorqueTime, periTorque,'.', 'Color',  [.5 .5 .5])
    line( [0 0], [-10 10]);
    xlim([-1 2]);
    ylim([-5 10]);
    title(Flickax, 'Overlayed Torque Triggers');
     
    % EMG Tile
    EMGax = nexttile(3);
    plot(EMGax,timeTQ, filtSolEMG,  'Color', unitColor(7,:))
    hold on
    plot(EMGax,timeTQ, filtMGEMG,  'Color', unitColor(5,:))
    plot(EMGax,smoothTime, smoothCST, 'Color',  [.5 .5 .5], 'LineWidth', 2)
    plot(EMGax,timeTQ, realEMG, 'Color', unitColor(1,:) ,'LineWidth', 2)
    xlim([0 length(smoothTime)/fsamp]);
    ylim([0 20]);
    legend('Sol', 'MG','TA CST', 'TA', 'location', 'NorthWest')
    ylabel('EMG Amplitude');
    title(EMGax, 'TA & Sol EMG & TA CST');
    
    % EMG Flicks Tile
    EMGflickax = nexttile(4);
    loopEMG = [];
    for i = 1:length(Flicks)
        loopMuscleEMGs = [];
        loopMuscleEMGs(1,:) = realEMG(timeTQ>Flicks(i)-2 & timeTQ<Flicks(i)+2);
        loopMuscleEMGs(2,:) = filtSolEMG(timeTQ>Flicks(i)-2 & timeTQ<Flicks(i)+2);
        loopEMG = [loopMuscleEMGs loopEMG];
    end
    
    EMG_Antagonist_Stim = loopEMG(2,:);
    EMG_Agonist_Stim = loopEMG(1,:);
    hold on
    plot(EMGflickax, tq_Stim_Time, EMG_Antagonist_Stim,'.','Color',  unitColor(7,:))
    plot(EMGflickax, tq_Stim_Time, EMG_Agonist_Stim,'.', 'Color', unitColor(1,:))
    line( [0 0], [-10 20]);
    xlim([-1 2]);
    ylim([0 20]);
    title(EMGflickax, 'Overlayed TA and Sol EMG Triggers');
    
    % PSTH Tile
    PSTHax = nexttile(5);
    bar(PSTHax, CST_Timebase,(CST_PSTH./(sum(CST_PSTH)/2)).*100, 'FaceColor', [.7 .7 .7], 'EdgeColor', [.7 .7 .7])
    hold on
    bar(PSTHax, Unit_Timebase,(Unit_PSTH./sum(Unit_PSTH)).*100,'FaceColor', unitColor(1,:),'EdgeColor', unitColor(1,:))
    xlim([-1 2]);
    yyaxis right ;
    plot(PSTHax,Unit_Timebase,cumsum(Unit_PSTH-mean(Unit_PSTH(Unit_Timebase<-0.5)))./length(Unit_PSTH),'LineWidth', 2)
    
    % PSF Tile
    PSFax = nexttile(6);
    title(PSFax, ' Discharge Rate and CST');
    hold on
    scatter(PSFax, CST_Stim_Time, CST_Stim_DR, 'filled', 'MarkerFaceColor', [.7 .7 .7])% '.', 'Color',  [.7 .7 .7], 'LineWidth', 10)
    scatter(Unit_Stim_Time(2:end), Unit_Stim_IDR, 'filled','MarkerFaceColor', unitColor(j,:), 'MarkerEdgeColor', 'k')
    line( [0 0], [-2 30]);
    ylim([-2 30]);
    xlim([-1 2]);
    ylabel('Discharge Rate');
    yyaxis right;
    cuSumTimebase = linspace(-2,2,length(Unit_Stim_IDR));
    plot(cuSumTimebase,cumsum(Unit_Stim_IDR-mean(Unit_Stim_IDR(Unit_Stim_Time(2:end)<-0.5)))./length(Unit_Stim_IDR),'LineWidth', 2)
    grid on;
    
    set(gcf, 'PaperPosition', [0 0 20 15]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [20 15]); %Set the paper to have width 5 and height 5.
    PSFFigname_save = fullfile(PSF_folder,[filename(1:end-30),'_Unit_',  num2str(j),'_PSF_PSTH']);
    print (gcf,'-dpdf',[PSFFigname_save '.pdf']);
    close(gcf);
    
end 

    
   