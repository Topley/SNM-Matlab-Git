function Pooled_Coherence_Analysis_v3_MT(fullFilename,Muscle1,Muscle2,Muscle3,Coco,fsamp,LW,LimFreq)
 close all %-except m1Fig m2Fig poolFig;
%  test file/variables for batch
%  filename = 'CoCotest02_60009_TA_1p_v23decomposed_MUCLEANED.mat';
%  Muscle1 = 'TA';
%  Muscle2 = 'Sol';
%  Muscle3 = 'MG';
%  N = 10;
%  fsamp = 2048;
%  LW = 1;
%  LimFreq = 500;
%  CoCo = 1;
% coherFolder = 'F:\Toolkit\Mirror\Coherence\DF Trials\Torque Feedback';
[folder, filename] = fileparts(fullFilename);
[~,coherLoc] = fileparts(folder);

if contains(folder, 'DF Trials')
    coherFolder = ['F:\Toolkit\Mirror\Coherence\DF Trials\' coherLoc];
else
    coherFolder = ['F:\Toolkit\Mirror\Coherence\PF Trials\' coherLoc];
end

if contains(folder, 'Coco')
    Coco = 2;   % 1 = single muscle, 2 = 2 muscles, 3 = 3 muscles;
else
    Coco = 1;
end

temp1 = load(fullFilename, 'MUPulses');
m1Pulses = temp1.MUPulses;

%%% sort units
m1MotorUnits = SortUnits(m1Pulses);

%Create agonist binary spike train for analysis and figures
[m1ST] = binarySpikeTrain(m1MotorUnits, []);
m1STFigure = m1ST;

% Ask if user wants to skip this file, prompts user to select analysis
% window by clicks
checkFig = figure;
clf('reset')
plot([0:length(m1STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m1STFigure'));
title(filename);
skip = input('skip trial? ','s')
if skip == 'y';
    return
end
disp('Choose analysis window ')
[ax1, ~] = ginput(2);
ax1 = round(ax1);

% If a multiple muscles are being analyzed, load other muscles file and MUs
% select contraction window based on range where motor units from all files
% are active
if Coco >= 2
    try
        m2Filename = replace(filename, Muscle1, Muscle2);
        temp2 = load(fullfile(folder,m2Filename), 'MUPulses');
        m2Pulses = temp2.MUPulses;
        m2MotorUnits = SortUnits(m2Pulses);
    catch
        disp('Could not find Muscle 2 file or Motor Units');
        try
            m2Filename = replace(filename, Muscle1, Muscle3);
        temp2 = load(fullfile(folder,m2Filename), 'MUPulses');
        m2Pulses = temp2.MUPulses;
        m2MotorUnits = SortUnits(m2Pulses);
        catch
            disp('Could not find second Muscle 2 file or Motor Units');
        end 
    end
    
    try
        [m2ST] = binarySpikeTrain(m2MotorUnits, []);
        m2STFigure = m2ST;
        [ax2] = analysisWindowClick(m2Filename,checkFig, m2STFigure);
        
        if ax1(1)>ax2(1)
            stax = ax1(1);
        else
            stax = ax2(1);
        end
        
        if ax1(2) < ax2(2)
            endax = ax1(2);
        else
            endax = ax2(2);
        end
        stax = round(stax);
        endax =round(endax);
        
    catch
        stax = round(ax1(1));
        endax = round(ax1(2));
        disp('Error plotting muscle 2 MUs');
    end
else
    stax = round(ax1(1));
    endax = round(ax1(2));
end

if Coco == 3
    try
        m3Filename = replace(filename, Muscle1, Muscle3);
        temp3 = load(m3Filename, 'MUPulses');
        m3Pulses = temp3.MUPulses;
        m3MotorUnits = SortUnits(m3Pulses);
    catch
        disp('Could not find Muscle 3 file or Motor Units');
    end
    
    try
        [m3ST] = binarySpikeTrain(m3MotorUnits, []);
        m3STFigure = m3ST;
        [ax3] = analysisWindowClick(m3Filename, m3STFigure);
        
        %Adjust to best range for both trials
        if ax3(1)> stax
            stax = ax3(1);
        end
        
        if ax3(2) < endax
            endax = ax3(2);
        end
        stax = round(stax);
        endax =round(endax);
        
    catch
        stax = round(ax1(1));
        endax = round(ax1(2));
        disp('Error plotting muscle 2 MUs');
    end
end

%%%% Getting proper axis range
axFig = figure;
a=axis;
a = [stax, endax, a(3), a(4)];
a = round(a*fsamp);
close(axFig)

%plot primary muscles spike trains to look for bad MUs & calculate COV
[m1goodUnits,m1goodFiring, m1COV, rSpike1] = removeBUsCoherence(filename, m1MotorUnits,m1ST,[stax endax]);
m1STFigure(rSpike1,:) = [];

%%%% Plot other muscles spike trains and remove bad ones & calculate COV
if Coco >= 2
    try
        [m2goodUnits,m2goodFiring, m2COV, rSpike2] = removeBUsCoherence(m2Filename, m2MotorUnits,m2ST,[stax endax]);
        m2STFigure(rSpike2,:) = [];
    catch
        disp('Cannot remove units from muscle #2 motor units');
    end
end

if Coco == 3
    try
        [m3goodUnits,m3goodFiring, m3COV, rSpike3] = removeBUsCoherence(m3Filename, m3MotorUnits,m3ST,[stax endax]);
        m3STFigure(rSpike3,:) = [];
    catch
        disp('Cannot remove units from muscle #3 motor units');
    end
end

try
%%%%  Replot smoothed STs without bad MUs and the average spike train
%     figure(3);
%     clf('reset')
% goodUnitFig = figure;
gcf;
 clf('reset')
hold on
    fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none');
%     hold on
    plot([0:length(m1STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m1STFigure'));
    plot([0:length(m1STFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');
%     
    if Coco >= 2
        try
            plot([0:length(m2STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m2STFigure'));
            plot([0:length(m2STFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),m2STFigure'),2),'Linewidth',3,'Color','k');
        catch
            disp('Cannot plot muscle #2 good motor units');
        end
    end
%     
%     if Coco == 3
%         try
%             plot(goodUnitFig,[0:length(m3STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m3STFigure'));
%             plot(goodUnitFig,[0:length(m3STFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),m3STFigure'),2),'Linewidth',3,'Color','k');
%         catch
%             disp('Cannot plot muscle #3 good motor units');
%         end
%     end
%     
catch
end
hold off


if Coco >= 2
%%%%% The first function (3 muscles) is incomplete and needs to be
    % finished, script will work fine for 1-2 muscles
    %       Suffix guide below
    %           r = random permutation of CSTs - intermuscular coherence
    %           p = all unique pairs of single MUs - intermuscular coherence
    %           ag = all unique pairs of agonist MUs  - intramuscular coherence
    %           ant = all unique pairs of antagonist MUs - intramuscular coherence
    try
        [cstCOHT, poolCOHT, m1COHT,m2COHT, cstF, poolF, m1F, m2F] = pooledCoherence3Muscle(m1goodFiring, m2goodFiring, m3goodFiring,LW,fsamp);
        
        pool3Firing = [m1goodFiring;m2goodFiring;m3goodFiring];
        m2COF = 1- (1-0.95)^(1/(size(afiring,2)/fsamp-1));
        poolCOF = 1- (1-0.95)^(1/(size(pool3Firing,2)/fsamp-1));
        m1COF = 1- (1-0.95)^(1/(size(firing,2)/fsamp-1));
        cstCOF = 1- (1-0.95)^(1/(size(firing,2)/fsamp-1));
        
        poolCOF2 = max(poolCOHT(poolF>LimFreq));
        m1COF2 = max(m1COHT(m1F>LimFreq));
        m2COF2 = max(m2COHT(m2F>LimFreq));
        cstCOF2 = max(cstCOHT(cstF>LimFreq));
        
        pooL = size(pool3Firing,2)/(fsamp*LW);
        m1L = size(firing,2)/(fsamp*LW);
        m2L = size(afiring,2)/(fsamp*LW);
        cstL = size(firing,2)/(fsamp*LW);
        
        poolZ = atanh(sqrt(poolCOHT))/sqrt(0.5/pooL);
        m1Z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L);
        m2Z = atanh(sqrt(m2COHT))/sqrt(0.5/m2L);
        cstZ = atanh(sqrt(cstCOHT))/sqrt(0.5/cstL);

    catch
        
%%%%%% Run function for intra and inter-muscular pooled coherence
        %       Suffix guide below
        %           cst = the sum of both muscles CSTs - intermuscular coherence
        %           p = all unique pairs of single MUs from boht muscles - intermuscular coherence
        %           ag = all unique pairs of agonist MUs  - intramuscular coherence
        %           ant = all unique pairs of antagonist MUs - intramuscular coherence
        
%         [poolFig,m1Fig,m2Fig,cstCOHT, poolCOHT, m1COHT,m2COHT, cstF, poolF, m1F, m2F] = pooledCoherence2MuscleTEST(m1goodFiring,m2goodFiring, LW,fsamp);
        [poolFig,m1Fig,m2Fig,cstCOHT, poolCOHT, m1COHT,m2COHT, cstF, poolF, m1F, m2F] = pooledCoherence2Muscle(m1goodFiring,m2goodFiring, LW,fsamp);
        keyboard
        poolFigname_save = fullfile(coherFolder,[filename(1:end-26) '_TAMG_Coherence_', num2str(stax), '_',num2str(endax)]);
            print (poolFig,'-dpdf',[poolFigname_save '.pdf']);
        m1Figname_save = fullfile(coherFolder,[filename(1:end-26) '_TATA_Coherence_', num2str(stax), '_',num2str(endax)]);
            print (m1Fig,'-dpdf',[m1Figname_save '.pdf']);
        m2Figname_save = fullfile(coherFolder,[filename(1:end-26) '_MGMG_Cohererence_', num2str(stax), '_',num2str(endax)]);
            print (m2Fig,'-dpdf',[m2Figname_save '.pdf']);
        
        pool2Firing = [m1goodFiring;m2goodFiring];
        m2COF = 1- (1-0.95)^(1/(size(m2goodFiring,2)/fsamp-1));
        poolCOF = 1- (1-0.95)^(1/(size(pool2Firing,2)/fsamp-1));
        m1COF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/fsamp-1));
        cstCOF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/fsamp-1));
        
        poolCOF2 = max(poolCOHT(poolF>LimFreq));
        m1COF2 = max(m1COHT(m1F>LimFreq));
        m2COF2 = max(m2COHT(m2F>LimFreq));
        cstCOF2 = max(cstCOHT(cstF>LimFreq));
        
        pooL = size(pool2Firing,2)/(fsamp*LW);
        m1L = size(m1goodFiring,2)/(fsamp*LW);
        m2L = size(m2goodFiring,2)/(fsamp*LW);
        cstL = size(m1goodFiring,2)/(fsamp*LW);
        
        poolZ = atanh(sqrt(poolCOHT))/sqrt(0.5/pooL);
        m1Z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L);
        m2Z = atanh(sqrt(m2COHT))/sqrt(0.5/m2L);
        cstZ = atanh(sqrt(cstCOHT))/sqrt(0.5/cstL);
        
    end
    
else
%%%%% Run pooled coherence for single muscle
    % Suffix guide below
    % r = random permutation of CSTs - intramuscular coherence
    % ag = all unique pairs of agonist MUs  - intramuscular coherence
    
    [m1Fig,m1COHT, cstCOHT, m1F, cstF] = pooledCoherence(m1goodFiring,LW,fsamp);
    %keyboard 
    m1Figname_save = fullfile(coherFolder,[filename(1:end-26) '_TATA_Coherence_', num2str(stax), '_',num2str(endax)]);
            print (m1Fig,'-dpdf',[m1Figname_save '.pdf']);
            
    m1COF = 1-(1-0.95)^(1/(size(m1goodFiring,2)/fsamp-1));
    m1COF2 = max(m1COHT(m1F>LimFreq));
    m1L = size(m1goodFiring,2)/(fsamp*LW);
    m1Z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L); % this is correct
    cstCOF = m1COF;
    cstCOF2 = max(cstCOHT(cstF>LimFreq));
    cstL = size(m1goodFiring,2)/(fsamp*LW);
    cstZ = atanh(sqrt(cstCOHT))/sqrt(0.5/cstL); % this is correct
    
end

% Summary top plot figure of smoothed STs with analysis window highlighted
% Antagonist units included if this is a coco trial
% removed STs grayed out
figure(101); clf('reset');subplot(2, 2, [1:2]);
hold all
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none')
plot([0:length(m1STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m1STFigure'));
plot([0:length(m1STFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');

if Coco >= 2
    try
        plot([0:length(m2STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m2STFigure'));
        plot([0:length(m2STFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),m2STFigure'),2),'Linewidth',3,'Color','k');
    end
end

if Coco == 3
    try
        plot([0:length(m3STFigure)-1]/fsamp,m3STFigure');
        plot([0:length(m3STFigure)-1]/fsamp,mean(m3STFigure',2),'Linewidth',3,'Color','k');
    end
end

hold off
title([filename '_' num2str(round(a(1)./fsamp)) '_' num2str(round(a(2)./fsamp))],'interpreter','none')
ylim([0 30])
xlabel('Time (s)')
ylabel('DR (pps)')

%%%% Summary coherence figure with the random perm and unique pairs
%%%% plotted. The 3 muscle plots need ot be added. This try section runs
%%%% 2 muscles

    try
        figure(101);
        subplot(2, 2, 3);
        plot(m1F,m1Z, 'SeriesIndex', 1)
        hold all
        plot(m2F,m2Z, 'SeriesIndex', 3)
        plot(poolF,poolZ, 'SeriesIndex', 5)
        m1COFz = atanh(sqrt(m1COF))/sqrt(0.5/m1L);
        h3 = plot(m1F,ones(1,length(m1F))*m1COFz,'r');
        m1COF2z = max(m1Z(m1F>LimFreq));
        h4 = plot(m1F,ones(1,length(m1F))*m1COF2z,'k','Linewidth',0.5);
        %hold off
        xlim([0 50])
        ylim([0 15])
        xlabel('Frequency (Hz)')
        ylabel('Coherence (z-transform)')
        m2COFz = atanh(sqrt(m2COF))/sqrt(0.5/m2L);
        m2COF2z = max(m2Z(m2F>LimFreq));
        poolCOFz = atanh(sqrt(poolCOF))/sqrt(0.5/pooL);
        poolCOF2z = max(m2Z(m2F>LimFreq));
%%%%% Second summary figure comparing the coherence plots of
        % Agonist intramuscular coherence
        % Antagonist intramuscular coherence
        % Agonist-Antagonist intermuscular coherence
        figure(99)
        t1 = tiledlayout(1,3);
        title(t1, 'All Unique Pairs Coherence')
        t1.XLabel.String = 'Frequency (Hz)';
        t1.YLabel.String = 'Coherence (z-transform)';
        nexttile(1)
        %     plot(rF,rz)
        plot(poolF,poolZ)
        hold on
        cstCOFz = atanh(sqrt(cstCOF))/sqrt(0.5/cstL);
        pCOFz = atanh(sqrt(poolCOF))/sqrt(0.5/pooL);
        plot(cstF,ones(1,length(cstF))*cstCOFz,'r');
        cstCOF2z = max(cstZ(cstF>LimFreq));
        pCOF2z = max(poolZ(poolF>LimFreq));
        plot(poolF,ones(1,length(poolF))*pCOF2z,'k','Linewidth',0.5);
        hold off
        title([Muscle1 '-' Muscle2])
        xlim([0 50])
        ylim([0 10])
        
        nexttile(2)
        plot(m1F,m1Z)
        hold on
        m1COFz = atanh(sqrt(m1COF))/sqrt(0.5/m1L);
        m1COF2z = max(m1Z(m1F>LimFreq));
        plot(m1F,ones(1,length(m1F))*m1COFz,'r');
        plot(m1F,ones(1,length(m1F))*m1COF2z,'k','Linewidth',0.5);
        hold off
        title([ Muscle1 '-' Muscle1])
        xlim([0 50])
        ylim([0 10])
        
        nexttile(3)
        plot(m2F,m2Z)
        hold on
        antCOFz = atanh(sqrt(m2COF))/sqrt(0.5/m2L);
        plot(m2F,ones(1,length(m2F))*antCOFz,'r');
        antCOF2z = max(m2Z(m2F>LimFreq));
        plot(m2F,ones(1,length(m2F))*antCOF2z,'k','Linewidth',0.5);
        hold off
        title([Muscle2 '-' Muscle2])
        xlim([0 50])
        ylim([0 10])
        
    catch
        
%%%% Summary coherence figure with the random perm and unique pairs
%%%% plotted. The 3 muscle plots need ot be added. This catch section runs
%%%% 1 muscle
        figure(101);
        subplot(2, 2, 3);
        plot(m1F,m1Z, 'SeriesIndex',1)
        hold on
        %     plot(rF, rz)
        m1COFz = atanh(sqrt(m1COF))/sqrt(0.5/m1L);
        cstCOFz = atanh(sqrt(cstCOF))/sqrt(0.5/cstL);
        plot(m1F,ones(1,length(m1F))*m1COFz,'r');
        m1COF2z = max(m1Z(m1F>LimFreq));
        cstCOF2z = max(cstZ(cstF>LimFreq));
        plot(m1F,ones(1,length(m1F))*m1COF2z,'k','Linewidth',0.5);
        hold off
        xlim([0 100])
        ylim([0 15])
        xlabel('Frequency (Hz)')
        ylabel('Coherence (z-transform)')
        
        % Second summary figure comparing the coherence plots of
        % Agonist intramuscular coherence
        % Antagonist intramuscular coherence
        % Agonist-Antagonist intermuscular coherence
        figure(99)
        t2 = tiledlayout(1,2);
        title(t2, 'CST vs pooled Intramuscular')
        t2.XLabel.String = 'Frequency (Hz)';
        t2.YLabel.String = 'Coherence (z-transform)';
        
        nexttile(1)
        plot(cstF,cstZ, 'k')
        hold on
        plot(cstF,ones(1,length(cstF))*cstCOFz,'r');
        plot(cstF,ones(1,length(cstF))*cstCOF2z,'k','Linewidth',0.5);
        hold off
        title('Random Perm CSTs')
        xlim([0 100])
        ylim([0 10])
        
        nexttile(2)
        plot(m1F,m1Z, 'Color', [.02 .02 .8])
        hold on
        m1COFz = atanh(sqrt(m1COF))/sqrt(0.5/m1L);
        m1COF2z = max(m1Z(m1F>LimFreq));
        plot(m1F,ones(1,length(m1F))*m1COFz,'r');
        plot(m1F,ones(1,length(m1F))*m1COF2z,'k','Linewidth',0.5);
        hold off
        title('All Unique Pairs')
        xlim([0 100])
        ylim([0 10])
        
    end

figure(101)
subplot(2, 2, 4); hold all
axis([0 40 0 40])

if a(2)>size(m1ST,2)
    a(2) =size(m1ST,2);
end

% Summary figure comparing meanDR vs CoVISI
% Agonist MUs Only
try
    subfiring = m1ST;
catch
    subfiring = m1ST(:,stax*fsamp:endax*fsamp);
end

for i = 1:size(subfiring,1)
    loopISI = diff(find(subfiring(i,:) == 1))./fsamp;
    loopISI = loopISI(loopISI>0.002 & loopISI<0.5);
    CoVISI(i) = std(loopISI)./mean(loopISI).*100;
    meanDR(i) = mean((1./loopISI));
    scatter(meanDR(i),CoVISI(i),'filled')
end

plottext = sprintf('meanDR = %.1f +/- %.0fpps \nCoVISI = %.0f +/- %.0f%',mean(meanDR),std(meanDR),mean(CoVISI),std(CoVISI));
text(1,3,plottext)
xlabel('meanDR (pps)')
ylabel('CoVISI (%)')


% List of variables saved
WinLen = LW*fsamp; % window length
%COHT % Coherence function
%F % frequencies
%fsamp % sampling frequency
%Iter % number of iterations
MUNumber = size(m1goodFiring,1); % number of motor units
MULength = size(m1goodFiring,2); % length of spike train in samp
%COF % standard confidence interval
%COF2 % peak COHT >500
%z % z transformed COHT
%COFz % standard z transform confidence interval
%COF2z  % peak z transform COHT >500
% a(1) % end
% a(2) % start
%rspike
% UsedMU = 1:(MUNumber+size(rspike))
% UsedMU(rspike) = []

if Coco >= 2
   try
        filename_save = [fullFilename(1:end-4) '_3muscle_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'m1COHT','m1F','m1COF','m1COF2','m1Z','m1COFz','m1COF2z',...
            'm2COHT','m2F','m2COF','m2COF2','m2Z','m2COFz','m2COF2z',...
            'm3COHT','m3F','m3COF','m3COF2','m3Z','m3COFz','m3COF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        
    catch
         filename_save = [fullFilename(1:end-4) '_2muscle_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'m1COHT','m1F','m1COF','m1COF2','m1Z','m1COFz','m1COF2z',...
            'm2COHT','m2F','m2COF','m2COF2','m2Z','m2COFz','m2COF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
    end
else
        filename_save = [fullFilename(1:end-4) '_1muscle_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'m1COHT','m1F','m1COF','m1COF2','m1Z','m1COFz','m1COF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);  
end 

end

function [m1Fig,m1COHT, cstCOHT, m1F, cstF] = pooledCoherence(firing,LW,fsamp)

% This function will take the smoothed binary spike trains of the active
% "good" units from ONE muscle and run the pooled coherence. 
% The autospectra of each unit is calculated, and so is the coherence 
% spectra of the two separate units. The autospectra of each unit is then
% subtracted out of the pooled coherence

N = 10;
%%%%%% Random permutations of active units
NREAL = size(firing,1);
group1 = [1:2:NREAL]; % take odds
group2 = [2:2:NREAL]; % take evens
if N<round(size(firing,1)/2-1)
    NI = N;
else
    NI = round(size(firing,1)/2-1);
end
h = waitbar(0)
PxxCST = 0;
PyyCST = 0;
PxyCST = 0;
Iter = 100;
for t = 1:Iter
    waitbar(t/Iter, h, 'Random Perm CST Coherence');
    [Pxx,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group1(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pyy,cstF] = cpsd(detrend(sum(firing(group2(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pxy,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    PxxCST = PxxCST + Pxx;
    PyyCST = PyyCST + Pyy;
    PxyCST = PxyCST + Pxy;
end
cstCOHT = abs(PxyCST).^2./(PxxCST.*PyyCST);

%%%%% All pairs of agonist coherence
m1Fig = figure('Visible', 'off');
m1Fax = axes('Parent',m1Fig);
hold(m1Fax,'on');
PxxM1 = 0;
PyyM1 = 0;
PxyM1 = 0;
unitL = [];
unitz = [];
unitCOHT =[];
loopCOHT =0;
% Iter = 100;
count = 0;
for j = 1:size(firing,1)-1
     waitbar(j/size(firing,1), h, 'Single Muscle All Pairs Coherence');

    for ij = j+1:size(firing,1)

    [Pxx,m1F] = cpsd(detrend(firing(j,:),0),detrend(firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pyy,m1F] = cpsd(detrend(firing(ij,:),0),detrend(firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pxy,m1F] = cpsd(detrend(firing(j,:),0),detrend(firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    PxxM1 = PxxM1 + Pxx;
    PyyM1 = PyyM1 + Pyy;
    PxyM1 = PxyM1 + Pxy;
    unitCOHT = abs(Pxy).^2./(Pxx.*Pyy);
    loopCOHT = loopCOHT+unitCOHT;
    %unitL = size(unitCOHT,2)/(fsamp*LW);
    unitL = 1-0.05^(1/(size(unitCOHT,1)/(fsamp*LW)-1));
    %unitz = atanh(sqrt(unitCOHT))/sqrt(0.5/unitL); % this is correct
    unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
    plot(m1Fax,m1F,unitz)
    count = count+1;
    end
end
close(h);
%m1COHT = loopCOHT./count;
m1COHT = abs(PxyM1).^2./(PxxM1.*PyyM1);
%testCOHT = mean(loopCOHT,2);
%m1L = size(m1COHT,2)/(fsamp*LW);
m1L = 1-0.05^(1/(size(m1COHT,1)/(fsamp*LW)-1));
%m1z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L); % this is correct
m1z = atanh(sqrt(m1COHT))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW))); % this is correct
bias = mean(m1z(m1F>250&m1F<500));
m1z = m1z-bias;
plot(m1Fax, m1F,m1z, 'Color','k', 'LineWidth', 2)

COF = 1-(1-0.95)^(1/(size(firing,2)/fsamp-1));
%COFz = atanh(sqrt(COF))/sqrt(0.5/m1L);
COFz = atanh(sqrt(COF))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW)))-bias; % this is correct
noiseCOF = max(m1COHT(m1F>500));
noiseCOFz = max(m1z(m1F>500));

plot(m1Fax,m1F,ones(1,length(m1F))*COFz,'r', 'LineWidth', 2);
plot(m1Fax,m1F,ones(1,length(m1F))*noiseCOFz,'b', 'LineWidth', 2);
xlim(m1Fax,[0 50]);
ylim(m1Fax,[0 1]);
yticks(m1Fax,[0 0.5 1]);
xlabel(m1Fax,'Frequency (Hz)')
ylabel(m1Fax,'Coherence (z-transform)')
title(m1Fax,'Pooled Coherence')
hold(m1Fax,'off');
set(m1Fig,'Visible','on');
%keyboard
end

function [poolFig,m1Fig,m2Fig,cstCOHT, pCOHT, m1COHT, m2COHT, cstF, pF, m1F, m2F] = pooledCoherence2Muscle(muscle1Firing,muscle2Firing,LW,fsamp)
% This function will take the smoothed binary spike trains of the active
% "good" units from TWO muscles and run the pooled coherence. 
% The autospectra of each unit is calculated from both muscles, and so is the coherence 
% spectra of the two separate units within and between muscles. 
% The autospectra of each unit is then subtracted out of the pooled
% coherence similar to the single muscle coherence

%%%% summed CST intermuscular coherence
h = waitbar(0);
try
PxxT1 = 0;
PyyT1 = 0;
PxyT1 = 0;
[Pxx,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle1Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[Pyy,cstF] = cpsd(detrend(sum(muscle2Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[Pxy,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);

PxxT1 = PxxT1 + Pxx;
PyyT1 = PyyT1 + Pyy;
PxyT1 = PxyT1 + Pxy;
cstCOHT = abs(PxyT1).^2./(PxxT1.*PyyT1);
catch
    disp('problem with CST coherence')
end 
%%%%%% 
% antagonist all permutations
%%%%%%
try
m2Fig = figure('Visible', 'off');
m2ax = axes('Parent',m2Fig);
hold(m2ax,'on');
Iterm2 = size(muscle2Firing,1);
Pxxm2 = 0;
Pyym2 = 0;
Pxym2 = 0;
unitL = [];
unitz = [];
unitCOHT = [];
loopCOHT =0;
count = 0;
for j = 1:size(muscle2Firing,1)-1
      waitbar(j/Iterm2,h, 'Processing Muscle 2 Coherence');
    for ij = j+1:size(muscle2Firing,1)
        [Pxx,m2F] = cpsd(detrend(muscle2Firing(j,:),0),detrend(muscle2Firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,m2F] = cpsd(detrend(muscle2Firing(ij,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,m2F] = cpsd(detrend(muscle2Firing(j,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        Pxxm2 = Pxxm2 + Pxx;
        Pyym2 = Pyym2 + Pyy;
        Pxym2 = Pxym2 + Pxy;
        unitCOHT = abs(Pxy).^2./(Pxx.*Pyy);
        loopCOHT = loopCOHT+unitCOHT;
        %unitL = size(unitCOHT,2)/(fsamp*LW);
        unitL = 1-0.05^(1/(size(unitCOHT,1)/(fsamp*LW)-1));
        %unitz = atanh(sqrt(unitCOHT))/sqrt(0.5/unitL); % this is correct
        unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
        plot(m2ax,m2F,unitz)
        count = count+1;
    end
end
m2COHT = abs(Pxym2).^2./(Pxxm2.*Pyym2);
%m2COHT = loopCOHT./count;
%m2L = 1-0.05^(1/(size(m2COHT,1)/(fsamp*LW)-1));
m2z = atanh(sqrt(m2COHT))/sqrt(0.5*(size(m2COHT,1)/(fsamp*LW))); % this is correct
biasM2 = mean(m2z(m2F>250 & m2F<500));
m2z = m2z - biasM2;
% m2L = size(m2COHT,2)/(fsamp*LW);
% m2z = atanh(sqrt(m2COHT))/sqrt(0.5/m2L); % this is correct
plot(m2ax, m2F,m2z, 'Color','k', 'LineWidth', 2)

m2COF = 1-(1-0.95)^(1/(size(muscle2Firing,2)/fsamp-1));
%m2COFz = atanh(sqrt(m2COF))/sqrt(0.5/m2L);
m2COFz = atanh(sqrt(m2COF))/sqrt(0.5*(size(m2COHT,1)/(fsamp*LW)))-biasM2; % this is correct
m2noiseCOF = max(m2COHT(m2F>500));
m2noiseCOFz = max(m2z(m2F>500));

plot(m2ax,m2F,ones(1,length(m2F))*m2COFz,'r', 'LineWidth', 2);
plot(m2ax,m2F,ones(1,length(m2F))*m2noiseCOFz,'b', 'LineWidth', 2);
xlim(m2ax,[0 40]);
ylim(m2ax,[0 1]);
yticks(m2ax,[0 .5 1]);
xlabel(m2ax,'Frequency (Hz)')
ylabel(m2ax,'Coherence (z-transform)')
title(m2ax,'MG-MG Coherence')
hold(m2ax,'off');

catch
    disp('problem with m2 coherence')
end 

%%%%%% 
% agonist all permutations
%%%%%%

try
m1Fig = figure('Visible', 'off');
m1Fax = axes('Parent',m1Fig);
hold(m1Fax,'on');
PxxM1 = 0;
PyyM1 = 0;
PxyM1 = 0;
unitL = [];
unitz = [];
loopCOHT = 0;
unitCOHT =[];
count =0;
Iterm1 = size(muscle1Firing,1);
for j = 1:size(muscle1Firing,1)-1
     waitbar(j/Iterm1,h, 'Processing Muscle 1 Coherence');
    for ij = j+1:size(muscle1Firing,1)
    [Pxx,m1F] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle1Firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pyy,m1F] = cpsd(detrend(muscle1Firing(ij,:),0),detrend(muscle1Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [Pxy,m1F] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle1Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    PxxM1 = PxxM1 + Pxx;
    PyyM1 = PyyM1 + Pyy;
    PxyM1 = PxyM1 + Pxy;
    unitCOHT = abs(Pxy).^2./(Pxx.*Pyy);
    loopCOHT = loopCOHT+unitCOHT;
     %unitL = size(unitCOHT,2)/(fsamp*LW);
    unitL = 1-0.05^(1/(size(unitCOHT,1)/(fsamp*LW)-1));
    %unitz = atanh(sqrt(unitCOHT))/sqrt(0.5/unitL); % this is correct
    unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
    plot(m1Fax,m1F,unitz)
    count = count+1;
    end
end
%m1COHT = loopCOHT./count;
%testCOHT = mean(loopCOHT,2);
%m1L = size(m1COHT,2)/(fsamp*LW);
m1COHT = abs(PxyM1).^2./(PxxM1.*PyyM1);
m1L = 1-0.05^(1/(size(m1COHT,1)/(fsamp*LW)-1));
%m1z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L); % this is correct
m1z = atanh(sqrt(m1COHT))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW))); % this is correct
bias = mean(m1z(m1F>250&m1F<500));
m1z = m1z-bias;
plot(m1Fax, m1F,m1z, 'Color','k', 'LineWidth', 2)

COF = 1-(1-0.95)^(1/(size(muscle1Firing,2)/fsamp-1));
%COFz = atanh(sqrt(COF))/sqrt(0.5/m1L);
COFz = atanh(sqrt(COF))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW)))-bias; % this is correct
noiseCOF = max(m1COHT(m1F>500));
noiseCOFz = max(m1z(m1F>500));

plot(m1Fax,m1F,ones(1,length(m1F))*COFz,'r', 'LineWidth', 2);
plot(m1Fax,m1F,ones(1,length(m1F))*noiseCOF,'b', 'LineWidth', 2);
xlim(m1Fax,[0 40]);
ylim(m1Fax,[0 1]);
yticks(m1Fax,[0 .5 1]);
xlabel(m1Fax,'Frequency (Hz)')
ylabel(m1Fax,'Coherence (z-transform)')
title(m1Fax,'TA-TA Coherence')
hold(m1Fax,'off');

catch
    disp('problem with m1 coherence')
end 
%%%%%%
% Pooled coherence between agonist and antagonist 
%%%%%%

try
poolFig = figure('Visible','off');
poolax = axes('Parent',poolFig);
hold(poolax,'on');
poolL = [];
poolz = [];

PxxT4 = 0;
PyyT4 = 0;
PxyT4 = 0;
unitCOHT = [];
loopCOHT =0;
count= 0;
Iterp = size(muscle1Firing,1);

for j = 1:size(muscle1Firing,1)
     waitbar(j/Iterp,h,'Processing Pooled Coherence');
    for ij = 1:size(muscle2Firing,1)
        [Pxx,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle1Firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,pF] = cpsd(detrend(muscle2Firing(ij,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT4 = PxxT4 + Pxx;
        PyyT4 = PyyT4 + Pyy;
        PxyT4 = PxyT4 + Pxy;
        unitCOHT = abs(Pxy).^2./((Pxx.*Pyy)+eps);
        loopCOHT = loopCOHT+unitCOHT;
        loopCOHT = loopCOHT+unitCOHT;
        pooL = 1-0.05^(1/(size(unitCOHT,1)/(fsamp*LW)-1));
        poolz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
        %poolL = size(unitCOHT,2)/(fsamp*LW);
        %poolz = atanh(sqrt(unitCOHT))/sqrt(0.5/poolL); % this is correct
        plot(poolax,pF,poolz)
        count = count+1;
    end
end
close(h);

pCOHT = abs(PxyT4).^2./(PxxT4.*PyyT4);
%pCOHT = loopCOHT./count;
pL = 1-0.05^(1/(size(pCOHT,1)/(fsamp*LW)-1));
pz = atanh(sqrt(pCOHT))/sqrt(0.5*(size(pCOHT,1)/(fsamp*LW))); % this is correct
biasPool = mean(pz(pF>250&pF<500));
pz = pz-biasPool;
%pL = size(pCOHT,2)/(fsamp*LW);
%pz = atanh(sqrt(pCOHT))/sqrt(0.5/pL); % this is correct
plot(poolax, pF,pz, 'Color','k', 'LineWidth', 2)

pCOF = 1-(1-0.95)^(1/(size(muscle1Firing,2)/fsamp-1));
%pCOFz = atanh(sqrt(pCOF))/sqrt(0.5/pL);
pCOFz = atanh(sqrt(pCOF))/sqrt(0.5*(size(pCOHT,1)/(fsamp*LW)))-biasPool; % this is correct
noiseIPCOF = max(pCOHT(pF>500));
noiseIPCOFz = max(pz(pF>500));

plot(poolax,pF,ones(1,length(pF))*pCOFz,'r', 'LineWidth', 2);
plot(poolax,pF,ones(1,length(pF))*noiseIPCOFz,'b', 'LineWidth', 2);
xlim(poolax,[0 40]);
ylim(poolax,[0 1]);
yticks(poolax,[0 .5 1]);
xlabel(poolax,'Frequency (Hz)')
ylabel(poolax,'Coherence (z-transform)')
title(poolax,'TA-MG Coherence')
hold(poolax,'off');

catch
    disp('problem with pooled coherence')
end 
set(m1Fig,'Visible','on');
set(m2Fig,'Visible','on');
set(poolFig,'Visible','on');

%pCOHT = abs(PxyT4).^2./((PxxT4.*PyyT4)+eps);
end

function [uniqueFig,m1Fig,m2Fig,cstCOHT, pCOHT, m1COHT, m2COHT, cstF, pcstF, m1F, m2F] = pooledCoherence2MuscleTEST(muscle1Firing,muscle2Firing,LW,fsamp)
% This function will take the smoothed binary spike trains of the active
% "good" units from TWO muscles and run the pooled coherence. 
% The autospectra of each unit is calculated from both muscles, and so is the coherence 
% spectra of the two separate units within and between muscles. 
% The autospectra of each unit is then subtracted out of the pooled
% coherence similar to the single muscle coherence

%%%% summed CST intermuscular coherence
try
concatenatedM1 = reshape(muscle1Firing,1,[]);
concatenatedM2 = reshape(muscle2Firing,1,[]);

m1CST = sum(muscle1Firing,1);
m1ReshapeCST =  repmat(m1CST,size(muscle2Firing,1));
concatencatedM1CST = m1ReshapeCST(1,:);

m2CST = sum(muscle2Firing,1);
m2ReshapeCST = repmat(m2CST,size(muscle1Firing,1));
concatenatedM2CST = m2ReshapeCST(1,:);
pooledCST = [concatencatedM1CST,concatenatedM2CST];

[PxxM1,m1F] = cpsd(detrend(concatenatedM1,0),detrend(concatenatedM1,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[PxxM2,m2F] = cpsd(detrend(concatenatedM2,0),detrend(concatenatedM2,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);

m1L = size(concatenatedM1,2)/(fsamp*LW);
m1z = atanh(sqrt(PxxM1))/sqrt(0.5/m1L); % this is correct
m2L = size(concatenatedM2,2)/(fsamp*LW);
m2z = atanh(sqrt(PxxM2))/sqrt(0.5/m2L); % this is correct

[PyyM1,m1CSTF] = cpsd(detrend(concatencatedM1CST,0),detrend(concatencatedM1CST,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[PyyM2,m2CSTF] = cpsd(detrend(concatenatedM2CST,0),detrend(concatenatedM2CST,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);

m1CSTL = size(concatencatedM1CST,2)/(fsamp*LW);
m1CSTz = atanh(sqrt(PyyM1))/sqrt(0.5/m1CSTL); % this is correct
m2CSTL = size(concatenatedM2CST,2)/(fsamp*LW);
m2CSTz = atanh(sqrt(PyyM2))/sqrt(0.5/m2CSTL); % this is correct

uniqueFig = figure('Visible','off');
uniqueax = axes('Parent',uniqueFig);
hold(uniqueax,'on');
plot(uniqueax,m1F,m1z)
plot(uniqueax,m1CSTF,m1CSTz)
plot(uniqueax,m2F,m2z)
plot(uniqueax,m2CSTF,m2CSTz)

COF = 1-(1-0.95)^(1/(size(muscle1Firing,2)/fsamp-1));
COFz = atanh(sqrt(COF))/sqrt(0.5/m1CSTL);
noiseCOF = max(m1CSTz(m1CSTF>500));
noiseCOFz = max(m1CSTz(m1CSTF>500));

plot(uniqueax,m1CSTF,ones(1,length(m1CSTF))*COFz,'r', 'LineWidth', 2);
plot(uniqueax,m1CSTF,ones(1,length(m1CSTF))*noiseCOFz,'b', 'LineWidth', 2);
xlim(uniqueax,[0 50]);
xlabel(uniqueax,'Frequency (Hz)')
ylabel(uniqueax,'Coherence (z-transform)')
title(uniqueax,'TATA & MGMG Total Coherence')
hold(uniqueax,'off');
set(uniqueFig,'Visible','on');

[PxyM1,uM1F] = cpsd(detrend(concatenatedM1,0),detrend(concatenatedM2CST,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[PxyM2,uM2F] = cpsd(detrend(concatenatedM2,0),detrend(concatencatedM1CST,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[Pxxyy,pcstF] = cpsd(detrend(pooledCST,0),detrend(pooledCST,0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);

uM1COHT = abs(PxyM1).^2./(PxxM1.*PyyM2);
uM2COHT = abs(PxyM2).^2./(PxxM2.*PyyM1);

m1UL = size(concatenatedM1,2)/(fsamp*LW);
m1Uz = atanh(sqrt(uM1COHT))/sqrt(0.5/m1UL); % this is correct
m2UL = size(concatenatedM2,2)/(fsamp*LW);
m2Uz = atanh(sqrt(uM2COHT))/sqrt(0.5/m2UL); % this is correct
pCSTL = size(pooledCST,2)/(fsamp*LW);
pCSTz = atanh(sqrt(Pxxyy))/sqrt(0.5/pCSTL); % this is correct

pCOF = 1-(1-0.95)^(1/(size(concatenatedM1,2)/fsamp-1));
pCOFz = atanh(sqrt(pCOF))/sqrt(0.5/m1UL);
pnoiseCOF = max(PxyM1(uM1F>500));
pnoiseCOFz = max(m1Uz(uM1F>500));

poolFig = figure('Visible','off');
poolax = axes('Parent',poolFig);
hold(poolax,'on');
%plot(poolax,pcstF,pCSTz)
plot(poolax,uM1F,m1Uz)
plot(poolax,uM2F,m2Uz)
plot(poolax,uM1F,ones(1,length(uM1F))*pCOFz,'r', 'LineWidth', 2);
plot(poolax,uM1F,ones(1,length(uM1F))*pnoiseCOFz,'b', 'LineWidth', 2);
xlim(poolax,[0 50]);
xlabel(poolax,'Frequency (Hz)')
ylabel(poolax,'Coherence (z-transform)')
title(poolax,'TATA & MGMG unique Coherence')
hold(poolax,'off');
set(poolFig,'Visible','on');

catch
    disp('problem with CST coherence')
end 

end

function [cstCOHT, poolCOHT, agCOHT, antCOHT, cstF, poolF, agF, antF] = pooledCoherence3Muscle(Agonsitfiring,LW,fsamp,Antagonistfiring)

%%%% summed CST intermuscular coherence
PxxT1 = 0;
PyyT1 = 0;
PxyT1 = 0;
[Pxx,cstF] = cpsd(detrend(sum(Agonsitfiring,1),0),detrend(sum(Agonsitfiring,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[Pyy,cstF] = cpsd(detrend(sum(Antagonistfiring,1),0),detrend(sum(Antagonistfiring,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
[Pxy,cstF] = cpsd(detrend(sum(Agonsitfiring,1),0),detrend(sum(Antagonistfiring,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);

PxxT1 = PxxT1 + Pxx;
PyyT1 = PyyT1 + Pyy;
PxyT1 = PxyT1 + Pxy;
cstCOHT = abs(PxyT1).^2./(PxxT1.*PyyT1);

%%%%%% muscle 1 all permutations

Iter = size(Antagonistfiring,1);
PxxT2 = 0;
PyyT2 = 0;
PxyT2 = 0;
for j = 1:size(Antagonistfiring,1)-1
    h= waitbar(j/Iter);
    for ij = j+1:size(Antagonistfiring,1)
        [Pxx,antF] = cpsd(detrend(Antagonistfiring(j,:),0),detrend(Antagonistfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,antF] = cpsd(detrend(Antagonistfiring(ij,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,antF] = cpsd(detrend(Antagonistfiring(j,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT2 = PxxT2 + Pxx;
        PyyT2 = PyyT2 + Pyy;
        PxyT2 = PxyT2 + Pxy;
    end
end
antCOHT = abs(PxyT2).^2./((PxxT2.*PyyT2)+eps);
close(h)

%%%%%% muscle 2 all permutations
Iter = size(Agonsitfiring,1);
PxxT3 = 0;
PyyT3 = 0;
PxyT3 = 0;
for j = 1:size(Agonsitfiring,1)-1
    h= waitbar(j/Iter);
    for ij = j+1:size(Agonsitfiring,1)
        [Pxx,agF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,agF] = cpsd(detrend(Agonsitfiring(ij,:),0),detrend(Agonsitfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,agF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT3 = PxxT3 + Pxx;
        PyyT3 = PyyT3 + Pyy;
        PxyT3 = PxyT3 + Pxy;
    end
end
agCOHT = abs(PxyT3).^2./((PxxT3.*PyyT3)+eps);
close(h)

%%%%%% muscle 3 all permutations
Iter = size(Agonsitfiring,1);
PxxT3 = 0;
PyyT3 = 0;
PxyT3 = 0;
for j = 1:size(Agonsitfiring,1)-1
    h= waitbar(j/Iter);
    for ij = j+1:size(Agonsitfiring,1)
        [Pxx,agF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,agF] = cpsd(detrend(Agonsitfiring(ij,:),0),detrend(Agonsitfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,agF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT3 = PxxT3 + Pxx;
        PyyT3 = PyyT3 + Pyy;
        PxyT3 = PxyT3 + Pxy;
    end
end
agCOHT = abs(PxyT3).^2./((PxxT3.*PyyT3)+eps);
close(h)

%%%% Pooled coherence between muscle 1 and muscle 2
Iter = size(Agonsitfiring,1);
PxxT4 = 0;
PyyT4 = 0;
PxyT4 = 0;
for j = 1:size(Agonsitfiring,1)
    h= waitbar(j/Iter);
    for ij = 1:size(Antagonistfiring,1)
        [Pxx,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,poolF] = cpsd(detrend(Antagonistfiring(ij,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT4 = PxxT4 + Pxx;
        PyyT4 = PyyT4 + Pyy;
        PxyT4 = PxyT4 + Pxy;
    end
end
poolCOHT = abs(PxyT4).^2./((PxxT4.*PyyT4)+eps);
close(h)

%%%% Pooled coherence between muscle 1 and muscle 3
Iter = size(Agonsitfiring,1);
PxxT4 = 0;
PyyT4 = 0;
PxyT4 = 0;
for j = 1:size(Agonsitfiring,1)
    h= waitbar(j/Iter);
    for ij = 1:size(Antagonistfiring,1)
        [Pxx,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,poolF] = cpsd(detrend(Antagonistfiring(ij,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT4 = PxxT4 + Pxx;
        PyyT4 = PyyT4 + Pyy;
        PxyT4 = PxyT4 + Pxy;
    end
end
poolCOHT = abs(PxyT4).^2./((PxxT4.*PyyT4)+eps);
close(h)

%%%% Pooled coherence between muscle 2 and muscle 3
Iter = size(Agonsitfiring,1);
PxxT4 = 0;
PyyT4 = 0;
PxyT4 = 0;
for j = 1:size(Agonsitfiring,1)
    h= waitbar(j/Iter);
    for ij = 1:size(Antagonistfiring,1)
        [Pxx,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Agonsitfiring(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,poolF] = cpsd(detrend(Antagonistfiring(ij,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,poolF] = cpsd(detrend(Agonsitfiring(j,:),0),detrend(Antagonistfiring(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxT4 = PxxT4 + Pxx;
        PyyT4 = PyyT4 + Pyy;
        PxyT4 = PxyT4 + Pxy;
    end
end
poolCOHT = abs(PxyT4).^2./((PxxT4.*PyyT4)+eps);
close(h)
end

function [] = savechoerence(firing,LW,fsamp,afiring)

end
