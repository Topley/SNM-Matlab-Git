function  Pooled_Coherence_Analysis_Base_LowerLeg(fullFilename,fsamp,LW,LimFreq)
 close all %-except m1Fig m2Fig poolFig;
 
%%%%  test file & variables for troubleshooting %%%%

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

temp1 = load(fullFilename, 'MUPulses');

% May need to load toraue to find the contraction hold time period
    % temp2 = load(fullFilename, 'Torque'); 
    % Torque = temp2.Torque;
    % m1Pulses = temp1.MUPulses;

%%% sort units
m1MotorUnits = SortUnits(m1Pulses);

%Create agonist binary spike train for analysis and figures
[m1ST] = binarySpikeTrain(m1MotorUnits, []);
m1STFigure = m1ST;

stax = 20;
endax = 50;

% If a multiple muscles are being analyzed, load other muscles file and MUs
% select contraction window based on range where motor units from all files
% are active

%%%% Getting proper axis range
axFig = figure;
a=axis;
a = [stax, endax, a(3), a(4)];
a = round(a*fsamp);
close(axFig)

%plot primary muscles spike trains to look for bad MUs & calculate COV
[m1goodUnits,m1goodFiring, m1COV, rSpike1] = removeBUsCoherence(filename, m1MotorUnits,m1ST,[stax endax]);
m1STFigure(rSpike1,:) = [];

try
%%%%  Replot smoothed STs without bad MUs and the average spike train
gcf;
 clf('reset')
hold on
    fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none');
    plot([0:length(m1STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m1STFigure'));
    plot([0:length(m1STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');
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
    
    [m1Fig, z, cstZ, F, cstF,freqband] = pooledCoherence(Muscle1, m1goodFiring,LW,fsamp);
    
    m1Figname_save = fullfile(coherFolder,[filename(1:end-26),'_', Muscle1,Muscle1,'_Coherence_', num2str(stax), '_',num2str(endax)]);
    print (m1Fig,'-dpdf',[m1Figname_save '.pdf']);
    
    COF = 1- (1-0.95)^(1/(size(m1goodFiring,2)/(fsamp*LW)-1));
    NoiseCL = max(z(F>LimFreq));
catch
    disp('Error in Muscle 1 coherence')
end


% Summary top plot figure of smoothed STs with analysis window highlighted
% removed STs grayed out
figure(101); clf('reset');subplot(2, 2, [1:2]);
hold all
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none')
plot([0:length(m1STFigure)-1]/fsamp,fftfilt(hanning(fsamp),m1STFigure'));
plot([0:length(m1STFigure)-1]/fsamp,mean(fftfilt(hanning(fsamp),m1STFigure'),2),'Linewidth',3,'Color','k');
hold off
title([filename '_' num2str(round(a(1)./fsamp)) '_' num2str(round(a(2)./fsamp))],'interpreter','none')
ylim([0 30])
xlabel('Time (s)')
ylabel('DR (pps)')

%%%% Summary coherence figure 

figure(101);
subplot(2, 2, 3);
plot(F,z, 'SeriesIndex', 1)
hold all
plot(F,ones(1,length(F))*COF,'r');
plot(F,ones(1,length(F))*NoiseCL,'k','Linewidth',0.5);
hold off
xlim([0 40])
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')

figure(101)
subplot(2, 2, 4); hold all
axis([0 40 0 40])

if a(2)>size(m1ST,2)
    a(2) =size(m1ST,2);
end

% Summary figure comparing meanDR vs CoVISI
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
    data.F = F;
    data.COF = COF;
    data.COF2 = NoiseCL;
    data.z = z;
    data.cstF = cstF;
    data.cstZ = cstZ;
    data.deltaband = freqband(1);
    data.alphaband = freqband(2);
    data.betaband = freqband(3);
    
    filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
    save([filename_save '.mat'],'data','fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
    print (gcf,'-dpdf',[filename_save '.pdf']);
    print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
end

end

function [m1Fig,m1z, cstCOHT, m1F, cstF, freqband] = pooledCoherence(muscle1, firing,LW,fsamp)

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
    unitz = [];
    unitCOHT = [];
    loopCOHT = 0;
    freqband = [];
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
            unitz = atanh(sqrt(unitCOHT))/sqrt(0.5*(size(unitCOHT,1)/(fsamp*LW)));
            plot(m1Fax,m1F,unitz)
        end
    end
    
    close(h);
    m1COHT = abs(PxyM1).^2./(PxxM1.*PyyM1);
    m1CL = 1-0.05^(1/(size(firing,2)/(fsamp*LW)-1));
    m1z = atanh(sqrt(m1COHT))/sqrt(0.5*(size(m1COHT,1)/(fsamp*LW))); % this is correct
    bias = mean(m1z(m1F>250&m1F<500));
    noiseCOFz = max(m1z(m1F>500));
    
    freqband(1) = trapz(m1F(1:50), m1z(1:50));    % 0 - 5Hz
    freqband(2) = trapz(m1F(50:150), m1z(50:150));    % 5 - 15Hz
    freqband(3) = trapz(m1F(150:350), m1z(150:350));    % 15 - 35Hz
    
    plot(m1Fax, m1F,m1z, 'Color','k', 'LineWidth', 2)
    area(m1Fax,m1F(1:50), m1z(1:50), 'SeriesIndex', 1, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(m1Fax,m1F(50:150), m1z(50:150), 'SeriesIndex', 3,'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(m1Fax,m1F(150:350), m1z(150:350),'SeriesIndex', 5, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
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
    
catch
    warning('Error processing single muscle all unique pairs intra-muscular coherence')
end
end

