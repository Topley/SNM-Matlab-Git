PSF_folder = 'F:\Toolkit\Mirror\CoCo09B\Decomposed\PSFs';

[filename,PathName] = uigetfile({'*.mat'},'Select the MATLAB Decomposed file');
load(filename);

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


TF = islocalmax(traceMax, 'MinSeparation', 5*fsamp, 'MinProminence', 2);
trigger = timeTQ(TF);
trigger(trigger<5) = [];
trigger(trigger>210) = [];
trigger = trigger-1;
realTorque = (TorqueFeedback-5).*2+6.5;
[CST_Stim_Time, CST_Stim_DR, CST_PSTH, CST_Timebase] = CST_PSF_PSTH(MUFiring, trigger)

for j = 1:length(MUFiring)
    
    avgDRs = [];
    unitFig = figure(j);
    set(unitFig,'Visible','off')
    t = tiledlayout(4,2);
    title(t,['Unit ', num2str(j)]) 
   
    
    ax1 = nexttile(1);
    hold on
    plot(timeTQ, realTorque, 'k')
    
    plot(timeTQ, traceReal, 'r')
    xlim([0 length(timeTQ)/fsamp])
    ylim([-10 20])
    ylabel('Torque')
    title(ax1, 'Torque Traces');
    
    ax3 = nexttile(3);
    plot(timeTQ, realEMG)
    hold on
    plot(smoothTime, smoothCST)
    xlim([0 length(smoothTime)/fsamp])
    ylim([0 20])
    ylabel('EMG Amplitude')
    title(ax3, 'TA EMG & CST');
    
    ax6 = nexttile(6);
    plot( CST_Stim_Time, CST_Stim_DR, '.', 'Color', 'k')
    
    for i = 1:length(trigger)
        mutime = MUFiring{j}./fsamp;
        muWindow = mutime>=trigger(i)-2.5 & mutime<=trigger(i)+2.5;
        DRs = mutime(muWindow);
        
        smoothWindow = smoothTime>=trigger(i)-2& smoothTime<=trigger(i)+2;
        eventWin = smoothTime(smoothWindow);
        figDRs = smoothIDRs(j,smoothWindow);
        
        torqwin = timeTQ>=trigger(i)-2 & timeTQ<=trigger(i)+2;
        torqWindow = timeTQ(torqwin);
        
        ax2 = nexttile(2);
        hold on
        grid on
        normTQTime = torqWindow-(trigger(i));
        normRealTime = torqWindow-(trigger(i)-.5);
        plot(ax2,normTQTime-.25,traceReal((trigger(i)-2)*fsamp:(trigger(i)+2)*fsamp),'r')
        plot(ax2,normRealTime-.25,realTorque((trigger(i)-2)*fsamp:(trigger(i)+2)*fsamp),'k')
        title(ax2, 'PF Triggers');
        line( [0 0], [-10 30])
        xlim([-1 2])
        ylim([-10 20])
        
        ax4 = nexttile(4);
        hold on
        plot(ax4,normTQTime,realEMG((trigger(i)-2)*fsamp:(trigger(i)+2)*fsamp),'k')
        grid on
        title(ax4, 'TA EMG Feedback');
        line( [0 0], [-5 30])
        xlim([-1 2])
        ylim([0 20])
        
        ax6 = nexttile(6);
        hold on
        normDRtime = DRs-(trigger(i));
        plot(ax6,normDRtime(2:end)+.25,1./diff(DRs), '.', 'MarkerSize', 12) ;
        %plot(ax7,normTQTime+.25,smoothCST((trigger(i)-2)*fsamp:(trigger(i)+2)*fsamp), 'Color', [.5 .5 .5])
        avgDRs = [avgDRs;figDRs];
        line( [0 0], [-2 30])
        ylim([-2 30])
        xlim([-1 2])
        ylabel('Discharge Rate')
        title(ax6, ' Discharge Rate and CST');
        grid on
        
        ax8 = nexttile(8);
        backGroundFiring = mutime>=trigger(i)-2.5 & mutime<=trigger(i);
        preFiring = mean(1./diff(mutime(backGroundFiring)));
        hold on
        test = 1./diff(DRs)-preFiring;
        grid on
        plot(normDRtime(2:end),cumsum(test)./preFiring)
        line( [0 0], [-10 10])
        xlim([-1 2])
        ylim([-10 10])
        xlabel('Peri-Stimulus Time')
        title(ax8, 'Cumulative Sum');
    end
    %plot(ax7,normTQTime,mean(avgDRs), 'LineWidth', 2)
    
    [psth, trialspx] = mpsth(MUFiring{j}./fsamp,trigger,'fr', 1, 'tb', 1, 'pre',1002,'post', 1002, 'binsz',2,'chart',0);
    nexttile(5)
    hold on
    bar(psth(:,1)+2,psth(:,2),'k','BarWidth',1)
    axis([min(psth(:,1))-10 max(psth(:,1))+10 0 max(psth(:,2))+1])
    ylabel(['counts per ' num2str(2) 'ms bin / fr (Hz)'])
    
    %     PSTHFigname_save = fullfile(PSF_folder,[filename(1:end-30),'_Unit_',  num2str(jj),'_PSTH']);
    
    nb_pulses_prior_trigger = sum(psth(find(psth<0),2)); % # of pulses occurring prior to triggering events
    baseline = nb_pulses_prior_trigger/length(find(psth(:,1)<0));
    cumulative_sum = cumsum(psth(:,2)-baseline); % baseline corrected cumulative sum
    nexttile(7)
    hold on
    plot(psth(:,1), cumulative_sum)
    xlim([-1000 1000])
    xlabel('Peri-Stimulus Time')
    ylabel('Normalized Cumulative Sum')
    set(unitFig,'Visible','on')
     
    set(gcf, 'PaperPosition', [0 0 20 15]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [20 15]); %Set the paper to have width 5 and height 5.
    PSFFigname_save = fullfile(PSF_folder,[filename(1:end-30),'_Unit_',  num2str(j),'_PSF_PSTH']);
    print (gcf,'-dpdf',[PSFFigname_save '.pdf']);
    close(gcf)
    
end