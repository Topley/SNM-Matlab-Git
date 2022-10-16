function Pooled_Coherence_Analysis_v4_MT(fullFilename, fsamp,LW,LimFreq, usermaxTorque)
close all %-except m1Fig m2Fig poolFig;
trial = 'Iramp';
    
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
[fileDir,filename,clusterFile, pdfdir] = setupDirectories(fullFilename, 'hold');

Muscle1 = 'TA';
Muscle2 = 'Sol';
Muscle3 = 'MG';

fileVariables = {'TraceFeedback', 'MUPulses', 'JR3Mz', 'EMG'};
load(fullfile(fileDir,filename), fileVariables{:});%, 'MUPulses');
displayPulses = MUPulses;
TraceFeedback = TraceFeedback;
EMGFeedback = filterEMG(EMG);

%%% Find the contraction hold window based on the EMG activity of muscle
%%% used for visual feedback
if trial == 'Iramp'
    checkFeedback = JR3Mz-mean(JR3Mz(1:fsamp));
    if mean(checkFeedback) < 2
        checkFeedback = (checkFeedback.*-1)./100;
        coco = 1;
    else
        checkFeedback = checkFeedback./100;
        coco = 0;
    end
    [stax] = findHolds(checkFeedback,EMGFeedback, coco);
    timeTQ = [0:length(JR3Mz)-1]/fsamp;
else
    checkFeedback = TraceFeedback-mean(TraceFeedback(1:fsamp));
    if mean(checkFeedback) < 2
        coco = 1;
    else
        coco = 0;
    end
    [stax] = findHolds(TraceFeedback,EMGFeedback, coco);
    timeTQ = [0:length(TraceFeedback)-1]/fsamp;
    
end

if contains(filename, 'Coco01')
    endax = stax+30;
else
    stax = stax+5;
    endax = stax+60;
end

[~, timeTorque, ~, cutTorqueFeedback, pfTQfeedback, dfTQfeedback] = TrimTorque([stax endax],JR3Mz, JR3Mz, usermaxTorque, coco);
torqueCOV = std(cutTorqueFeedback./mean(cutTorqueFeedback))*100;

%%% sort units
displayMotorUnits = SortUnits(displayPulses);

%Create agonist binary spike train for analysis and figures
[displayST] = binarySpikeTrain(displayMotorUnits, []);
displaySTFigure = displayST;

%%% Plot units with window for user to check
checkFig = figure;
clf('reset')
hold on
fill([stax, endax,endax,stax],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none');
plot([0:length(displaySTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),displaySTFigure'));
title(filename);
hold off

% If a multiple muscles are being analyzed, load other muscles file and MUs
% select contraction window based on range where motor units from all files
% are active

[TAUnits,SolUnits, MGUnits,TAEMGFeedback,SolEMGFeedback,MGEMGFeedback] = getAllCoherenceUnits(fullfile(fileDir,filename),clusterFile, coco);
muscleList = [];

if size(TAUnits,2) >4
    try
        [TAST, ~] = binarySpikeTrain(TAUnits, JR3Mz);
        TASTFigure = TAST;
        [TAgoodUnits,spiketrainTA,TAMUFiring, TACOV, TArSpike] = removeBUsCoherenceAuto(TAUnits,TAST,[stax endax]);
        TASTFigure(TArSpike,:) = [];
        muscleList = [muscleList; 'TA'];
        torqueCorrTA = corrcoef(sum(fftfilt(hanning(fsamp),TAMUFiring'),2), cutTorqueFeedback(fsamp:end));
    end
end

if size(MGUnits,2) >4
    try
        [MGST] = binarySpikeTrain(MGUnits, JR3Mz);
        MGSTFigure = MGST;
        [MGgoodUnits,MGMUFiring, MGCOV, MGrSpike] = removeBUsCoherenceAuto(MGUnits,MGST,[stax endax]);
        MGSTFigure(MGrSpike,:) = [];
        muscleList = [muscleList; 'MG'];
        torqueCorrMG = corrcoef(sum(fftfilt(hanning(fsamp),MGMUFiring'),2), cutTorqueFeedback(fsamp:end));
    end
end

if size(SolUnits,2) >4
    try
        [SolST] = binarySpikeTrain(SolUnits, JR3Mz);
        SolSTFigure = SolST;
        [SolgoodUnits,SolMUFiring, SolCOV, SolrSpike] = removeBUsCoherenceAuto(SolUnits,SolST,[stax endax]);
        SolSTFigure(SolrSpike,:) = [];
        muscleList = [muscleList; 'Sol'];
        torqueCorrSol = corrcoef(sum(fftfilt(hanning(fsamp),SolMUFiring'),2), cutTorqueFeedback(fsamp:end));
    end
end

try
    TAMGMUFiring = [TAMUFiring;MGMUFiring];
    torqueCorrTAMG = corrcoef(sum(fftfilt(hanning(fsamp),TAMGMUFiring'),2), cutTorqueFeedback(fsamp:end));
end 

try
    TASolMUFiring = [TAMUFiring;SolMUFiring];
    torqueCorrTASol = corrcoef(sum(fftfilt(hanning(fsamp),TASolMUFiring'),2), cutTorqueFeedback(fsamp:end));
end 

try
    SolMGMUFiring = [SolMUFiring;MGMUFiring];
    torqueCorrSolMG = corrcoef(sum(fftfilt(hanning(fsamp),SolMGMUFiring'),2), cutTorqueFeedback(fsamp:end));
end 


%%%% Getting proper axis range
axFig = figure;
a=axis;
a = [stax, endax, a(3), a(4)];
a = round(a*fsamp);
close(axFig)

trimSaveFilename = fullfile(fileDir, filename(1:12));
filename_save = [trimSaveFilename '_ALLUnits_pooled_coherence_' num2str(stax) '_' num2str(endax),'.mat'];
% if exist(filename_save)
%     return
% end

%%%%  Replot smoothed STs without bad MUs and the average spike train
%     figure(3);
%     clf('reset')
% goodUnitFig = figure;
cMap = Tableau;
checkFig;
clf('reset');
hold on
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none');
try
    plot([0:length(TASTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),TASTFigure'),'Color', cMap(2,:));
    plot([0:length(TASTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),TASTFigure'),2),'Linewidth',3,'Color', cMap(1,:));
end
try
    plot([0:length(SolSTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),SolSTFigure'),'Color', cMap(4,:));
    plot([0:length(SolSTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),SolSTFigure'),2),'Linewidth',3, 'Color', cMap(3,:));
end

try
    plot([0:length(MGSTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),MGSTFigure'),'Color', cMap(6,:));
    plot([0:length(MGSTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),MGSTFigure'),2),'Linewidth',3,'Color', cMap(5,:));
    
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
    [TAFig, TAZ, TAcstZ, TAF, TAcstF,TAfreqband,TAexplained] = pooledIntraCoherence(Muscle1, TAMUFiring,LW,fsamp);
    
    TAFigname_save = fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle1,'-',Muscle1,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (TAFig,'-dpdf',[TAFigname_save '.pdf']);
    
    TACOF = 1- (1-0.95)^(1/(size(TAMUFiring,2)/(fsamp*LW)-1));
    TANoiseCL = max(TAZ(TAF>LimFreq));

catch
    disp([Muscle1, 'coherence unable to be processed'])
end

try
    [SolFig, SolZ, SolcstZ, SolF, SolcstF,Solfreqband, Solexplained] = pooledIntraCoherence(Muscle2, SolMUFiring,LW,fsamp);
    
    SolFigname_save =fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle2,'-',Muscle2,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (SolFig,'-dpdf',[SolFigname_save '.pdf']);
    
    SolCOF = 1- (1-0.95)^(1/(size(SolMUFiring,2)/(fsamp*LW)-1));
    SolNoiseCL = max(SolZ(SolF>LimFreq));
    
catch
    disp([Muscle2, 'coherence unable to be processed'])
end

try
    [TASolpoolFig,TASolcstCOHT, TASolpoolZ,TASolcstF, TASolpoolF,TASolexplained] = pooledInterCoherence(Muscle1, Muscle2,TAMUFiring, SolMUFiring,LW,fsamp);
    
    TASolpoolFigname_save = fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle1,'-',Muscle2,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (TASolpoolFig,'-dpdf',[TASolpoolFigname_save '.pdf']);
    
    TASolPoolCOF = 1- (1-0.95)^(1/(size(TAMUFiring,2)/(fsamp*LW)-1));
    TASolPoolNoiseCL = max(TASolpoolZ(TASolpoolF>LimFreq));
    
catch
    disp([Muscle1, Muscle2, 'coherence unable to be processed'])
end

try
    [MGFig, MGZ, MGcstZ, MGF, MGcstF,MGexplained] = pooledIntraCoherence(Muscle3, MGMUFiring,LW,fsamp);
    
    MGFigname_save = fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle3,'-',Muscle3,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (MGFig,'-dpdf',[MGFigname_save '.pdf']);
    
    MGCOF = 1- (1-0.95)^(1/(size(MGMUFiring,2)/(fsamp*LW)-1));
    MGNoiseCL = max(MGZ(MGF>LimFreq));
    
catch
    disp([Muscle3, 'coherence unable to be processed'])
end

try
    [TAMGpoolFig,TAMGcstCOHT, TAMGpoolZ,TAMGcstF, TAMGpoolF,TAMGexplained] = pooledInterCoherence(Muscle1, Muscle3,TAMUFiring, MGMUFiring,LW,fsamp);
    
    TAMGpoolFigname_save = fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle1,'-',Muscle3,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (TAMGpoolFig,'-dpdf',[TAMGpoolFigname_save '.pdf']);
    
    TAMGPoolCOF = 1- (1-0.95)^(1/(size(TAMUFiring,2)/(fsamp*LW)-1));
    TAMGPoolNoiseCL = max(TAMGpoolZ(TAMGpoolF>LimFreq));
catch
    disp([Muscle1, Muscle3, 'coherence unable to be processed'])
end
try
    [SolMGpoolFig,SolMGcstCOHT, SolMGpoolZ,SolMGcstF, SolMGpoolF,SolMGexplained] = pooledInterCoherence(Muscle2, Muscle3,SolMUFiring, MGMUFiring,LW,fsamp);
    
    SolMGpoolFigname_save = fullfile(fileDir,'coherence_pdf\',[filename(1:12),'_', Muscle2,'-',Muscle3,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (SolMGpoolFig,'-dpdf',[SolMGpoolFigname_save '.pdf']);
    
    SolMGPoolCOF = 1- (1-0.95)^(1/(size(SolMUFiring,2)/(fsamp*LW)-1));
    SolMGPoolNoiseCL = max(SolMGpoolZ(SolMGpoolF>LimFreq));
catch
    disp([Muscle2, Muscle3, 'coherence unable to be processed'])
end

% Summary top plot figure of smoothed STs with analysis window highlighted
% Antagonist units included if this is a coco trial
% removed STs grayed out
TorqueFeedback = (JR3Mz-(mean(JR3Mz(1:fsamp))));
newTorqueFeedback = TorqueFeedback./100;

if mean(newTorqueFeedback) < 0 && coco == 0
    newTorqueFeedback = newTorqueFeedback.*-1;
end

figure(101); clf('reset');subplot(2, 2, [1:2]);
hold all
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none')
plot([0:length(newTorqueFeedback)-1]/fsamp,newTorqueFeedback, 'k')
try
    plot([0:length(TASTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),TASTFigure'),'Color', cMap(2,:));
end

try
    plot([0:length(SolSTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),SolSTFigure'),'Color', cMap(4,:));
end

try
    plot([0:length(MGSTFigure)-1]/fsamp,2*fftfilt(hanning(fsamp),MGSTFigure'),'Color', cMap(6,:));
end

try
    plot([0:length(TASTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),TASTFigure'),2),'Linewidth',3,'Color', cMap(1,:));
end

try
    plot([0:length(SolSTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),SolSTFigure'),2),'Linewidth',3, 'Color', cMap(3,:));
end

try
    plot([0:length(MGSTFigure)-1]/fsamp,mean(2*fftfilt(hanning(fsamp),MGSTFigure'),2),'Linewidth',3,'Color', cMap(5,:));
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
hold all
try
    plot(TAF,TAZ, 'SeriesIndex', 1)
end
try
    plot(SolF,SolZ, 'SeriesIndex', 3)
end
try
    plot(MGF,MGZ, 'SeriesIndex', 5)
end

try
    plot(TASolpoolF,TASolpoolZ, 'k')
end
try
    plot(TAMGpoolF,TAMGpoolZ, 'g')
end
try
    plot(SolMGpoolF,SolMGpoolZ, 'r')
end

try
    plot(TAF,ones(1,length(TAF))*TACOF,'r');
    plot(TAF,ones(1,length(TAF))*TANoiseCL,'k','Linewidth',0.5);
catch
    try
        plot(SolF,ones(1,length(SolF))*SolCOF,'r');
        plot(SolF,ones(1,length(SolF))*SolNoiseCL,'k','Linewidth',0.5);
    catch
        plot(MGF,ones(1,length(MGF))*MGCOF,'r');
        plot(MGF,ones(1,length(MGF))*MGNoiseCL,'k','Linewidth',0.5);
    end
end

hold off
xlim([0 40])
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')


%%%% Summary coherence figure with the random perm and unique pairs
%%%% plotted. The 3 muscle plots need ot be added. This catch section runs
%%%% 1 muscle

figure(101)
subplot(2, 2, 4); hold all
axis([0 40 0 40])

if a(2)>size(displayST,2)
    a(2) =size(displayST,2);
end

% Summary figure comparing meanDR vs CoVISI
% Agonist MUs Only
try
    TAsubfiring = TAMUFiring;
catch
    try
        TAsubfiring = TAMUFiring(:,stax*fsamp:endax*fsamp);
    catch
        TAsubfiring = [];
    end
end

try
    Solsubfiring = SolMUFiring;
catch
    try
        Solsubfiring = SolMUFiring(:,stax*fsamp:endax*fsamp);
    catch
        Solsubfiring = [];
    end
end

try
    MGsubfiring = MGMUFiring;
catch
    try
        MGsubfiring = MGMUFiring(:,stax*fsamp:endax*fsamp);
    catch
        MGsubfiring = [];
    end
end

try
    for i = 1:size(TAsubfiring,1)
        loopISI = diff(find(TAsubfiring(i,:) == 1))./fsamp;
        loopISI = loopISI(loopISI>0.002 & loopISI<0.5);
        TACoVISI(i) = std(loopISI)./mean(loopISI).*100;
        TAmeanDR(i) = mean((1./loopISI));
        scatter(TAmeanDR(i),TACoVISI(i),'filled', 'SeriesIndex', 1)
    end
end

hold on
try
    for i = 1:size(Solsubfiring,1)
        loopISI = diff(find(Solsubfiring(i,:) == 1))./fsamp;
        loopISI = loopISI(loopISI>0.002 & loopISI<0.5);
        SolCoVISI(i) = std(loopISI)./mean(loopISI).*100;
        SolmeanDR(i) = mean((1./loopISI));
        scatter(SolmeanDR(i),SolCoVISI(i),'filled','SeriesIndex', 3)
    end
end

try
    for i = 1:size(MGsubfiring,1)
        loopISI = diff(find(MGsubfiring(i,:) == 1))./fsamp;
        loopISI = loopISI(loopISI>0.002 & loopISI<0.5);
        MGCoVISI(i) = std(loopISI)./mean(loopISI).*100;
        MGmeanDR(i) = mean((1./loopISI));
        scatter(MGmeanDR(i),MGCoVISI(i),'filled','SeriesIndex', 3)
    end
end

if contains(fullFilename, 'TA')
    plottext = sprintf('meanDR = %.1f +/- %.0fpps \nCoVISI = %.0f +/- %.0f%',mean(TAmeanDR),std(TAmeanDR),mean(TACoVISI),std(TACoVISI));
    MUNumber = size(TAMUFiring,1); % number of motor units
    MULength = size(TAMUFiring,2); % length of spike train in samp
elseif contains(fullFilename, 'Sol')
    plottext = sprintf('meanDR = %.1f +/- %.0fpps \nCoVISI = %.0f +/- %.0f%',mean(SolmeanDR),std(SolmeanDR),mean(SolCoVISI),std(SolCoVISI));
    MUNumber = size(SolMUFiring,1); % number of motor units
    MULength = size(SolMUFiring,2); % length of spike train in samp
else
    plottext = sprintf('meanDR = %.1f +/- %.0fpps \nCoVISI = %.0f +/- %.0f%',mean(MGmeanDR),std(MGmeanDR),mean(MGCoVISI),std(MGCoVISI));
    MUNumber = size(MGMUFiring,1); % number of motor units
    MULength = size(MGMUFiring,2); % length of spike train in samp
end
text(1,4,plottext)
xlabel('meanDR (pps)')
ylabel('CoVISI (%)')


% List of variables saved
WinLen = LW*fsamp; % window length
%COHT % Coherence function
%F % frequencies
%fsamp % sampling frequency
%Iter % number of iterations
try
    data.TAMUNumber = size(TAMUFiring,1); % number of motor units
    data.TAMULength = size(TAMUFiring,2); % length of spike train in samp
end
try
    data.SolMUNumber = size(SolMUFiring,1); % number of motor units
    data.SolMULength = size(SolMUFiring,2); % length of spike train in samp
end
try
    data.MGMUNumber = size(MGMUFiring,1); % number of motor units
    data.MGMULength = size(MGMUFiring,2); % length of spike train in samp
end
% MUNumber = size(TAMUFiring,1); % number of motor units
% MULength = size(TAMUFiring,2); % length of spike train in samp
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

try
    data.TAF = TAF;
    data.TACOF = TACOF;
    data.TACOF2z = TANoiseCL;
    data.TAZ = TAZ;
    data.TAcstF = TAcstF;
    data.TAcstZ = TAcstZ;
    data.TAdeltaband = TAfreqband(1);
    data.TAalphaband = TAfreqband(2);
    data.TAbetaband = TAfreqband(3);
    data.TAexplained = TAexplained;
    data.torqueCorrTA = torqueCorrTA;
end

try
    data.TASolZ = TASolpoolZ;
    data.TASolF = TASolpoolF;
    data.TASolcstF = TASolcstF;
    data.TASolcstZ = TASolcstCOHT;
    data.TASolexplained = TASolexplained;
end

try
    data.TAMGZ = TAMGpoolZ;
    data.TAMGF = TAMGpoolF;
    data.TAMGcstF = TAMGcstF;
    data.TAMGcstZ = TAMGcstCOHT;
    data.TAMGexplained = TAMGexplained;
end


try
    data.MGF = MGF;
    data.MGZ = MGZ;
    data.MGCOF = MGCOF;
    data.MGCOF2z = MGNoiseCL;
    data.MGcstF = MGcstF;
    data.MGcstZ = MGcstZ;
    data.MGdeltaband = MGfreqband(1);
    data.MGalphaband = MGfreqband(2);
    data.MGbetaband = MGfreqband(3);
    data.MGexplained = MGexplained;
    data.torqueCorrMG = torqueCorrMG;
end

try
    data.SolF = SolF;
    data.SolZ = SolZ;
    data.SolCOF = SolCOF;
    data.SolCOF2z = SolNoiseCL;
    data.SolcstF = SolcstF;
    data.SolcstZ = SolcstZ;
    data.Soldeltaband = Solfreqband(1);
    data.Solalphaband = Solfreqband(2);
    data.Solbetaband = Solfreqband(3);
    data.torqueCorrSol = torqueCorrSol;
end

try
    data.SolMGZ = SolMGpoolZ;
    data.SolMGF = SolMGpoolF;
    data.SolMGcstF = SolMGcstF;
    data.SolMGcstZ = SolMGcstCOHT;
end


save(filename_save,'data','fsamp', 'WinLen','a','torqueCOV');
print (figure(101),'-dpdf',[filename_save(1:end-4) '.pdf']);
end

function [m1Fig,m1z, cstCOHT, m1F, cstF, freqband,explained] = pooledIntraCoherence(muscle1, firing,LW,fsamp)

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
catch
    disp('Error processing single muscle all unique pairs intra-muscular coherence')
end
close(h);
end

function [poolFig,cstCOHT, poolZ, cstF, pF,explained] = pooledInterCoherence(muscle1, muscle2, muscle1Firing,muscle2Firing,LW,fsamp)
% This function will take the smoothed binary spike trains of the active
% "good" units from TWO muscles and run the pooled coherence.
% The autospectra of each unit is calculated from both muscles, and so is the coherence
% spectra of the two separate units within and between muscles.
% The autospectra of each unit is then subtracted out of the pooled
% coherence similar to the single muscle coherence
% wcoherence(sum(muscle1Firing,1), sum(muscle2Firing,1), fsamp)
%%%% summed CST intermuscular coherence
firing = [muscle1Firing;muscle2Firing];
[coeff,score,latent,tsquared,explained] = pca(rot90(2*fftfilt(hanning(fsamp),firing')));
[lambda,psi,T,stats,F] = factoran(2*fftfilt(hanning(fsamp),firing'),1, 'rotate','orthomax');

[lambda,psi,T,stats,F] = factoran(smoothST,3, 'rotate','orthomax');

h = waitbar(0, ['Processing ', muscle1, '-',muscle2,' CST Coherence']);
try
    [M1xx,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle1Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [M2yy,cstF] = cpsd(detrend(sum(muscle2Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    [M12xy,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
    
    cstCOHT = abs(M12xy).^2./(M1xx.*M2yy);
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
% set(poolFig,'Visible','on');
close(h);
end