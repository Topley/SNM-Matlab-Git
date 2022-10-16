function Pooled_Coherence_Analysis_v4_MT(fullFilename,fsamp,LW,LimFreq)
 close all %-except m1Fig m2Fig poolFig;
 
%%%%  test file/variables for batch %%%%

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

%%%%%

[folder, filename] = fileparts(fullFilename);
[~,coherLoc] = fileparts(folder);

if contains(folder, 'DF Trials')
    coherFolder = ['F:\Toolkit\Mirror\Coherence\DF Trials\' coherLoc];
    Muscle1 = 'TA';
    Muscle2 = 'MG';
    Muscle3 = 'Sol';
    Muscle4 = 'LG';
else
    coherFolder = ['F:\Toolkit\Mirror\Coherence\PF Trials\' coherLoc];
    Muscle1 = 'MG';
    Muscle2 = 'TA';
    Muscle3 = 'Sol';
    Muscle4 = 'LG';
end

if contains(folder, 'Coco')
    Coco = 2;   % 1 = single muscle, 2 = 2 muscles, 3 = 3 muscles;
else
    Coco = 1;
end

temp1 = load(fullFilename, 'MUPulses');
m1Pulses = temp1.MUPulses;
temp2 = load(fullFilename, 'TorqueFeedback');
TorqueFeedback = movmean(temp2.TorqueFeedback,fsamp);
timeAxes = [0:length(TorqueFeedback)-1]./fsamp;
TorqueChanges = ischange(TorqueFeedback,'linear','MaxNumChanges',4);
timeChange = timeAxes(TorqueChanges);

stax = round(timeChange(2));
endax = round(timeChange(2));


%%% sort units
m1MotorUnits = SortUnits(m1Pulses);

%Create agonist binary spike train for analysis and figures
[m1ST] = binarySpikeTrain(m1MotorUnits, []);
m1STFigure = m1ST;

% Ask if user wants to skip this file, prompts user to select analysis
% window by clicks
% checkFig = figure;
% clf('reset')
% plot([0:length(m1STFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),m1STFigure'));
% title({filename; coherLoc});
% skip = input('skip trial? ','s')
% if skip == 'y';
%     return
% end
% disp('Choose analysis window ')
% [ax1, ~] = ginput(2);
% ax1 = round(ax1);

stax = ax1(1);
endax = ax1(2);

% If a multiple muscles are being analyzed, load other muscles file and MUs
% select contraction window based on range where motor units from all files
% are active

if contains(Muscle1,'TA')
    try
        m2Filename = replace(filename, Muscle1, Muscle2);
        temp2 = load(fullfile(folder,m2Filename), 'MUPulses');
        m2Pulses = temp2.MUPulses;
        m2MotorUnits = SortUnits(m2Pulses);
        Coco = 2;
    catch
        disp('No MG MUCLEANED file');
    end
    try
        m3Filename = replace(filename, Muscle1, Muscle3);
        temp3 = load(fullfile(folder,m3Filename), 'MUPulses');
        m3Pulses = temp3.MUPulses;
        m3MotorUnits = SortUnits(m3Pulses);
        Coco = 3;
    catch
        disp('No Sol MUCLEANED file');
    end
    
elseif contains(Muscle1,'MG')
    try
        m2Filename = replace(filename, Muscle1, Muscle2);
        temp2 = load(fullfile(folder,m2Filename), 'MUPulses');
        m2Pulses = temp2.MUPulses;
        m2MotorUnits = SortUnits(m2Pulses);
        Coco = 2;
    catch
        disp('No TA MUCLEANED file');
    end
    
    try
        m3Filename = replace(filename, Muscle1, Muscle3);
        temp3 = load(fullfile(folder,m3Filename), 'MUPulses');
        m3Pulses = temp3.MUPulses;
        m3MotorUnits = SortUnits(m3Pulses);
        Coco = 3;
    catch
        disp('No Sol MUCLEANED file');
    end
end

try
    [m2ST] = binarySpikeTrain(m2MotorUnits);
    m2STFigure = m2ST;
    [ax2] = analysisWindowClick(m2Filename,checkFig, m2STFigure,[stax endax]);
    ax2 = round(ax2);
    
    if stax<ax2(1)
        stax = ax2(1);
    end
    
    if endax > ax2(2)
        endax = ax2(2);
    end
    
catch
    disp('Error plotting muscle 2 MUs');
end


try
    [m3ST] = binarySpikeTrain(m3MotorUnits);
    m3STFigure = m3ST;
    [ax3] = analysisWindowClick(m3Filename, checkFig, m3STFigure,[stax endax]);
    ax3 = round(ax3);
    %Adjust to best range for both trials
    if ax3(1)> stax
        stax = ax3(1);
    end
    
    if  endax > ax3(2)
        endax = ax3(2);
    end
    
catch
    disp('Error plotting muscle 3 MUs');
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
    plot([0:length(m1STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m1STFigure'));
    plot([0:length(m1STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');
%     
    if Coco >= 2
        try
            plot([0:length(m2STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m2STFigure'));
            plot([0:length(m2STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m2STFigure'),2),'Linewidth',3,'Color','k');
        catch
            disp('Cannot plot muscle #2 good motor units');
        end
    end
    
    if Coco == 3
        try
            plot([0:length(m3STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m3STFigure'));
            plot([0:length(m3STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m3STFigure'),2),'Linewidth',3,'Color','k');
        catch
            disp('Cannot plot muscle #3 good motor units');
        end
    end
    
catch
end
hold off

%%%%% The first function (3 muscles) is incomplete and needs to be
    % finished, script will work fine for 1-2 muscles
    %       Suffix guide below
    %           r = random permutation of CSTs - intermuscular coherence
    %           p = all unique pairs of single MUs - intermuscular coherence
    %           ag = all unique pairs of agonist MUs  - intramuscular coherence
    %           ant = all unique pairs of antagonist MUs - intramuscular coherence

    try
        
        [m1Fig, m1Z, cstZ, m1F, cstF,freqband] = pooledIntraCoherence(Muscle1, m1goodFiring,LW,fsamp);
        
        m1Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle1,Muscle1,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m1Fig,'-dpdf',[m1Figname_save '.pdf']);
        
        m1COF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/(fsamp*LW)-1));
        m1NoiseCL = max(m1Z(m1F>LimFreq));
    catch
        disp('Error in Muscle 1 coherence')
    end
    try
        [m2Fig, m2Z, m2cstZ, m2F, m2cstF,freqband] = pooledIntraCoherence(Muscle2, m2goodFiring,LW,fsamp);

        m2Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle2,Muscle2,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m2Fig,'-dpdf',[m2Figname_save '.pdf']);
        
        m2COF = 1- (1-0.95)^(1/(size(m2goodFiring,2)/(fsamp*LW)-1));
        m2NoiseCL = max(m2Z(m2F>LimFreq));
    catch
        disp('Error in Muscle 2 coherence')
    end
    try
        [m3Fig, m3Z, m3cstZ, m3F, m3cstF] = pooledIntraCoherence(Muscle3, m3goodFiring,LW,fsamp);
        m3Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle3,Muscle3,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m3Fig,'-dpdf',[m3Figname_save '.pdf']);
        
        m3COF = 1- (1-0.95)^(1/(size(m3goodFiring,2)/(fsamp*LW)-1));
        m3NoiseCL = max(m3Z(m3F>LimFreq));
    catch
        disp('Error in Muscle 3 coherence')
    end
    
    try
        [m12poolFig,m12cstCOHT, m12poolZ,m12cstF, m12poolF,] = pooledInterCoherence(Muscle1, Muscle2,m1goodFiring, m2goodFiring,LW,fsamp);
        pool12Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle1,Muscle2,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m12poolFig,'-dpdf',[pool12Figname_save '.pdf']);
        m12PoolCOF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/(fsamp*LW)-1));
        m12PoolNoiseCL = max(m12poolZ(m12poolF>LimFreq));
    catch
        disp('Error in Muscle 1 & 2 pooled coherence')
    end
    
    try
         [m13poolFig,m13cstCOHT, m13poolZ,m13cstF, m13poolF,] = pooledInterCoherence(Muscle1, Muscle3,m1goodFiring, m3goodFiring,LW,fsamp);
        pool13Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle1,Muscle3,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m13poolFig,'-dpdf',[pool13Figname_save '.pdf']);
        m13PoolCOF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/(fsamp*LW)-1));
        m13PoolNoiseCL = max(m13poolZ(m13poolF>LimFreq));
    catch
        disp('Error in Muscle 1 & 3 pooled coherence')
    end
    try
        [m23poolFig,m23cstCOHT, m23poolZ,m23cstF, m23poolF,] = pooledInterCoherence(Muscle2, Muscle3,m2goodFiring, m3goodFiring,LW,fsamp);
        pool23Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle2,Muscle3,'_Coherence_', num2str(stax), '_',num2str(endax)]);
        print (m23poolFig,'-dpdf',[pool23Figname_save '.pdf']);
        m23PoolCOF = 1- (1-0.95)^(1/(size(m2goodFiring,2)/(fsamp*LW)-1));
        m23PoolNoiseCL = max(m23poolZ(m23poolF>LimFreq));
    catch
        disp('Error in Muscle 2 & 3 pooled coherence')
    end
%         m1COF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/(fsamp*LW)-1));
%         m1NoiseCL = max(poolZ(poolF>LimFreq));

% Summary top plot figure of smoothed STs with analysis window highlighted
% Antagonist units included if this is a coco trial
% removed STs grayed out
figure(101); clf('reset');subplot(2, 2, [1:2]);
hold all
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none')
plot([0:length(m1STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m1STFigure'));
plot([0:length(m1STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');

if Coco >= 2
    try
        plot([0:length(m2STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m2STFigure'));
        plot([0:length(m2STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m2STFigure'),2),'Linewidth',3,'Color','k');
    end
end

if Coco == 3
    try
       plot([0:length(m3STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m3STFigure'));
        plot([0:length(m3STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m3STFigure'),2),'Linewidth',3,'Color','k');
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

figure(101);
subplot(2, 2, 3);
plot(m1F,m1Z, 'SeriesIndex', 1)
hold all
try
    plot(m2F,m2Z, 'SeriesIndex', 3)
    plot(m12poolF,m12poolZ, 'k')
end
plot(m1F,ones(1,length(m1F))*m1COF,'r');
plot(m1F,ones(1,length(m1F))*m1NoiseCL,'k','Linewidth',0.5);
hold off
xlim([0 40])
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')
        
%%%%% Second summary figure comparing the coherence plots of
        % Agonist intramuscular coherence
        % Antagonist intramuscular coherence
        % Agonist-Antagonist intermuscular coherence
        try
            figure(99)
            t1 = tiledlayout(1,3);
            title(t1, 'All Unique Pairs Coherence')
            t1.XLabel.String = 'Frequency (Hz)';
            t1.YLabel.String = 'Coherence (z-transform)';
            nexttile(1)
            plot(m12poolF,m12poolZ)
            title([Muscle1 '-' Muscle2])
            xlim([0 40])
            ylim([0 1])
            
            nexttile(2)
            plot(m1F,m1Z)
            title([ Muscle1 '-' Muscle1])
            xlim([0 40])
            ylim([0 1])
            
            nexttile(3)
            plot(m2F,m2Z)
            title([Muscle2 '-' Muscle2])
            xlim([0 40])
            ylim([0 1])
            
        catch
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
            hold off
            title('Random Perm CSTs')
            xlim([0 40])
            ylim([0 1])
            
            nexttile(2)
            plot(m1F,m1Z, 'Color', [.02 .02 .8])
            hold on
            plot(m1F,ones(1,length(m1F))*m1COF,'r');
            plot(m1F,ones(1,length(m1F))*m1NoiseCL,'k','Linewidth',0.5);
            hold off
            title('All Unique Pairs')
            xlim([0 40])
            ylim([0 1])
        end
%%%% Summary coherence figure with the random perm and unique pairs
%%%% plotted. The 3 muscle plots need ot be added. This catch section runs
%%%% 1 muscle

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

if contains(filename, 'TA')
    data.TAF = m1F;
    data.TACOF = m1COF;
    data.TACOF2z = m1NoiseCL;
    data.TAZ = m1Z;
    data.cstF = cstF;
    data.cstZ = cstZ;
    data.deltaband = freqband(1);
    data.alphaband = freqband(2);
    data.betaband = freqband(3);
    
    try
        data.pool12Z = m12poolZ;
        data.pool12F = m12poolF;
        data.cst12F = m12cstF; 
        data.cst12Z = m12cstCOHT;
    end 
    
    try 
        data.MGF = m2F;
        data.MGZ = m2Z;
        data.MGCOF = m2COF;
        data.MGCOF2z = m2NoiseCL;
    end
    
    try
        data.SolF = m3F;
        data.SolZ = m3Z;
        data.SolCOF = m3COF;
        data.SolCOF2z = m3NoiseCL;
    end
    
    filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
    save([filename_save '.mat'],'data','fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
    print (gcf,'-dpdf',[filename_save '.pdf']);
    print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
    
elseif contains(filename, 'MG')
    data.MGF = m1F;
    data.MGCOF = m1COF;
    data.MGCOF2 = m1NoiseCL;
    data.MGZ = m1Z;
    
    try 
        data.TAF = m2F;
        data.TAZ = m2Z;
        data.TACOF = m2COF;
        data.TACOF2z = m2NoiseCL;
    end
    
    try    
        data.SolF = m3F;
        data.SolZ = m3Z;
        data.SolCOF = m3COF;
        data.SolCOF2z = m3NoiseCL;
    end
    
    filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
    save([filename_save '.mat'],'data','fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
    print (gcf,'-dpdf',[filename_save '.pdf']);
    print(figure(99),'-dpdf',[filename_save '_UnitPlots.pdf']);
    
else
    data.SolF = m1F;
    data.SolCOF = m1COF;
    data.SolCOF2 = m1NoiseCL;
    data.SolZ = m1Z;
    
    try 
        data.TAF = m2F;
        data.TAZ = m2Z;
         data.TACOF = m2COF;
        data.TACOF2z = m2NoiseCL;
    end
    
    try
        data.MGF = m3F;
        data.MGZ = m3Z;
        data.MGCOF = m3COF;
        data.MGCOF2z = m3NoiseCL;
    end
    
    filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
    save([filename_save '.mat'],'data','fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
    print (gcf,'-dpdf',[filename_save '.pdf']);
    print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
end

end

function [m1Fig,m1z, cstCOHT, m1F, cstF, freqband] = pooledIntraCoherence(muscle1, firing,LW,fsamp)

% This function will take the smoothed binary spike trains of the active
% "good" units from ONE muscle and run the pooled coherence. 
% The autospectra of each unit is calculated, and so is the coherence 
% spectra of the two separate units. The autospectra of each unit is then
% subtracted out of the pooled coherence
try
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
    
    h = waitbar(0, 'Processing 100-Perm  CST Coherence');
    PxxCST = 0;
    PyyCST = 0;
    PxyCST = 0;
    Iter = 100;
    for t = 1:Iter
        [Pxx,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group1(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,cstF] = cpsd(detrend(sum(firing(group2(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxCST = PxxCST + Pxx;
        PyyCST = PyyCST + Pyy;
        PxyCST = PxyCST + Pxy;
    end
    cstCOHT = abs(PxyCST).^2./(PxxCST.*PyyCST);
catch
    disp('Error processing single muscle CST intra-muscular coherence')
end

%%%%% All pairs of agonist coherence
try
    m1Fig = figure('Visible', 'off');
    m1Fax = axes('Parent',m1Fig);
    hold(m1Fax,'on');
    PxxM1 = 0;
    PyyM1 = 0;
    PxyM1 = 0;
    %unitL = [];
    unitz = [];
    unitCOHT =[];
    loopCOHT =0;
    % Iter = 100;
    waitbar(0.5, h, ['Processing ', muscle1, ' Units']);
    
    for j = 1:size(firing,1)-1
        
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
            %unitCL = 1-0.05^(1/(size(unitCOHT,1)/(fsamp*LW)-1));
            %unitz = atanh(sqrt(unitCOHT))/sqrt(0.5/unitL); % this is correct
            unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
            plot(m1Fax,m1F,unitz)
        end
    end
    
    close(h);
    % keyboard
    m1COHT = abs(PxyM1).^2./(PxxM1.*PyyM1);
    %m1L = size(m1COHT,2)/(fsamp*LW);
    %m1z = atanh(sqrt(m1COHT))/sqrt(0.5/m1L); % this is correct
    m1CL = 1-0.05^(1/(size(firing,2)/(fsamp*LW)-1));
    m1z = atanh(sqrt(m1COHT))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW))); % this is correct
    bias = mean(m1z(m1F>250&m1F<500));
    freqband = [];
    freqband(1) = trapz(m1F(1:50), m1z(1:50));    % 0 - 5Hz
    freqband(2) = trapz(m1F(50:150), m1z(50:150));    % 5 - 15Hz
    freqband(3) = trapz(m1F(150:350), m1z(150:350));    % 15 - 35Hz
    %m1z = m1z-bias;
    plot(m1Fax, m1F,m1z, 'Color','k', 'LineWidth', 2)
    area(m1Fax,m1F(1:50), m1z(1:50), 'SeriesIndex', 1, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(m1Fax,m1F(50:150), m1z(50:150), 'SeriesIndex', 3,'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(m1Fax,m1F(150:350), m1z(150:350),'SeriesIndex', 5, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    
    %m1CL = m1CL-bias;
    %COF = 1-(1-0.95)^(1/(size(firing,2)/fsamp-1));
    %COFz = atanh(sqrt(COF))/sqrt(0.5/m1L);
    %COFz = atanh(sqrt(COF))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW)))-bias; % this is correct
    %noiseCOF = max(m1COHT(m1F>500));
    noiseCOFz = max(m1z(m1F>500));
    
    plot(m1Fax,m1F,ones(1,length(m1F))*m1CL,'r', 'LineWidth', 2);
    plot(m1Fax,m1F,ones(1,length(m1F))*noiseCOFz,'b', 'LineWidth', 2);
    xlim(m1Fax,[0 40]);
    ylim(m1Fax,[0 1]);
    yticks(m1Fax,[0 0.5 1]);
    xlabel(m1Fax,'Frequency (Hz)')
    ylabel(m1Fax,'Coherence (z-transform)')
    title(m1Fax,[muscle1,'-', muscle1,' Pooled Coherence'])
    hold(m1Fax,'off');
    set(m1Fig,'Visible','on');
    %keyboard
catch
    disp('Error processing single muscle all unique pairs intra-muscular coherence')
end
end

function [poolFig,cstCOHT, poolZ, cstF, pF] = pooledInterCoherence2Muscle(muscle1, muscle2, muscle1Firing,muscle2Firing,LW,fsamp)
% This function will take the smoothed binary spike trains of the active
% "good" units from TWO muscles and run the pooled coherence. 
% The autospectra of each unit is calculated from both muscles, and so is the coherence 
% spectra of the two separate units within and between muscles. 
% The autospectra of each unit is then subtracted out of the pooled
% coherence similar to the single muscle coherence

%%%% summed CST intermuscular coherence
h = waitbar(0, ['Processing ', muscle1, '-',muscle2,' CST Coherence']);
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
% Pooled coherence between agonist and antagonist 
%%%%%%

try
poolFig = figure('Visible','off');
poolax = axes('Parent',poolFig);
hold(poolax,'on');
PxxPool = 0;
PyyPool = 0;
PxyPool = 0;
Unitz = [];
unitCOHT = [];
loopCOHT = 0;
Iterp = size(muscle1Firing,1);
waitbar(0.5,h,'Processing Pooled Coherence');
  
for j = 1:size(muscle1Firing,1)
    
    for ij = 1:size(muscle2Firing,1)
        [Pxx,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle1Firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pyy,pF] = cpsd(detrend(muscle2Firing(ij,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        [Pxy,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
        PxxPool = PxxPool + Pxx;
        PyyPool = PyyPool + Pyy;
        PxyPool = PxyPool + Pxy;
        unitCOHT = abs(Pxy).^2./((Pxx.*Pyy)+eps);
        loopCOHT = loopCOHT+unitCOHT;
        Unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW))); % this is correct
        plot(poolax,pF,Unitz)
    end
end

close(h);
pCOHT = abs(PxyPool).^2./(PxxPool.*PyyPool);
pCL = 1-0.05^(1/(size(muscle1Firing,2)/(fsamp*LW)-1));
poolZ = atanh(sqrt(pCOHT))/sqrt(0.5*(size(pCOHT,1)/(fsamp*LW))); % this is correct
biasPool = mean(poolZ(pF>250&pF<500));
poolZ = poolZ-biasPool;
pCL = pCL-biasPool;
noisePoolCOFz = max(poolZ(pF>500));

plot(poolax, pF,poolZ, 'Color','k', 'LineWidth', 2)
plot(poolax,pF,ones(1,length(pF))*pCL,'r', 'LineWidth', 2);
plot(poolax,pF,ones(1,length(pF))*noisePoolCOFz,'b', 'LineWidth', 2);

xlim(poolax,[0 40]);
ylim(poolax,[0 1]);
yticks(poolax,[0 .5 1]);
xlabel(poolax,'Frequency (Hz)')
ylabel(poolax,'Coherence (z-transform)')
title(poolax,[muscle1,'-',muscle2, ' Coherence'])
hold(poolax,'off');
catch
    disp('problem with pooled coherence')
end 
set(poolFig,'Visible','on');
pause(2)
end