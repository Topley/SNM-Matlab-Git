function Pooled_Coherence_Analysis_Base(fullFilename, fsamp,LW,LimFreq, usermaxTorque)
close all %-except m1Fig m2Fig poolFig;


%%%%  test file/variables for batch %%%%
%  filename = 'Coco01_13210_TA_5_v23decomposed_MUCLEANED.mat';
%  N = 10;
%  fsamp = 2048;
%  LW = 1;
%  LimFreq = 500;

%% Gets all files associated with trial and sets up folder to print pdf figures
[fileDir,Cleanedfile,Decompfile, Clusterfile, pdfdir] = Set_Directories(fullFilename, 'hold');

%% Extract units, trace, EMG, and torque
cleanVars = matfile(fullfile(fileDir,Cleanedfile));
clusterVars = matfile(fullfile(fileDir,Clusterfile));
MUPulses = cleanVars.MUPulses;
TargetTrace = clusterVars.TraceFeedback;

basedTrace = TargetTrace-mean(TargetTrace(1:fsamp));
TargetFeedback = movmean(basedTrace,fsamp)./10;

EMGFeedback = EMGpro(clusterVars.EMG, 'channel', 27, 'rmsWin', 500);

rawTorque = clusterVars.Torque;
rawTorqueFeedback = clusterVars.JR3Mz;
basedTorqueFeedback = (rawTorqueFeedback-mean(rawTorqueFeedback(1:fsamp)))./100;
basedTorque = (rawTorque-mean(rawTorque(1:fsamp)))./100;
if mean(basedTorqueFeedback) < 1
    basedTorqueFeedback = basedTorqueFeedback .*-1;
end

if mean(basedTorque) < 1
    basedTorque = basedTorque.*-1;
end

newTorqueFeedback = movmean(basedTorqueFeedback,fsamp);
newTorque = basedTorque;

%% Find the contraction hold window based on feedback or EMG activity of muscle used for visual feedback
stax = Set_HoldWindow(TargetFeedback,EMGFeedback);
endax = stax + 30;
timeTQ = [0:length(newTorqueFeedback)-1]/fsamp;

%% Trim torque to the hold window
cutTorqueFeedback = newTorqueFeedback((stax+1)*fsamp:endax*fsamp);
torqueCOV = std(cutTorqueFeedback./mean(cutTorqueFeedback))*100;

%% sort and remove units, generate binary spike train and smoothed spike train, get COV of smoothed units (for removal), ...
% get correlation between smooth spt and torque
MUPulses = SortUnits(MUPulses);
[binSPT, smoothBinSPT] = BinarySpikeTrain(MUPulses, newTorqueFeedback);
figSPT = smoothBinSPT;
[cutMUFiring, cutMU_SPT ,cutSmoothSPT, rSpikeCOV, rSpikes] = Remove_MUs_Auto(MUPulses,binSPT,[stax endax]);
figSPT(rSpikes,:) = [];
smoothCST = sum(cutSmoothSPT)./size(cutSmoothSPT,1);
torqueCorr = corrcoef(smoothCST, cutTorqueFeedback);

%% Getting proper axis range
axFig = figure;
a=axis;
a = [stax, endax, a(3), a(4)];
a = round(a*fsamp);
close(axFig)

trimSaveFilename = fullfile(fileDir, Cleanedfile(1:12));
filename_save = [trimSaveFilename '_ALLUnits_pooled_coherence_' num2str(stax) '_' num2str(endax),'.mat'];

cMap = Tableau;

newStr = extractBetween(filename,"_R","_1_v23");
Muscle1 = newStr;
try
    [cFig, m1Z, m1F, coherCI, NoiseCI, freqband] = Pooled_Intramuscular_Coherence(Muscle1, cutMU_SPT,1,2048);
    
    %blue line > noise over 500 Hz, Red line > than CI based off all unit
    %comparisons
    Figname_save = fullfile(fileDir,pdfdir,[Cleanedfile(1:12),'_', Muscle1,'-',Muscle1,'_coherence_', num2str(stax), '_',num2str(endax)]);
    print (cFig,'-dpdf',[Figname_save '.pdf']);
catch
    disp([Muscle1, 'coherence unable to be processed'])
end

% Summary top plot figure of smoothed STs with analysis window highlighted
% Antagonist units included if this is a coco trial
% removed STs grayed out
figure(101); clf('reset');subplot(2, 2, [1:2]);
hold all
fill([a(1)./fsamp, a(2)./fsamp,a(2)./fsamp,a(1)./fsamp],[0 0 30 30],'k','facealpha',0.1,'edgecolor','none')
plot([0:length(newTorqueFeedback)-1]/fsamp,newTorqueFeedback, 'k')

plot([0:length(figSPT)-1]/fsamp,figSPT','Color', cMap(2,:));
plot([0:length(figSPT)-1]/fsamp,mean(figSPT,2),'Linewidth',3,'Color', cMap(1,:));
hold off

title([Cleanedfile '_' num2str(round(a(1)./fsamp)) '_' num2str(round(a(2)./fsamp))],'interpreter','none')
ylim([0 30])
xlabel('Time (s)')
ylabel('DR (pps)')

%% Summary coherence figure with the random perm and unique pairs

figSPT(101);
subplot(2, 2, 3);
hold all

plot(m1F,m1Z, 'SeriesIndex', 1)
plot(m1F,ones(1,length(m1F))*coherCI,'r');
plot(m1F,ones(1,length(m1F))*NoiseCI,'k','Linewidth',0.5);

hold off
xlim([0 40])
ylim([0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')

figSPT(101)
subplot(2, 2, 4); hold all
axis([0 40 0 40])

if a(2)>size(displaySpT,2)
    a(2) =size(displaySpT,2);
end

% Summary figure comparing meanDR vs CoVISI
% Agonist MUs Only
try
    subfiring = cutMUFiring;
catch
    try
        subfiring = cutMUFiring(:,stax*fsamp:endax*fsamp);
    catch
        disp('Cant plot COVISI')
    end
end

hold on
for i = 1:size(subfiring,1)
    loopISI = diff(find(subfiring(i,:) == 1))./fsamp;
    loopISI = loopISI(loopISI>0.002 & loopISI<0.5);
    CoVISI(i) = std(loopISI)./mean(loopISI).*100;
    meanDR(i) = mean((1./loopISI));
    scatter(meanDR(i),CoVISI(i),'filled', 'SeriesIndex', 1)
end



plottext = sprintf('meanDR = %.1f +/- %.0fpps \nCoVISI = %.0f +/- %.0f%',mean(meanDR),std(meanDR),mean(CoVISI),std(CoVISI));
MUNumber = size(cutMUFiring,1); % number of motor units
MULength = size(cutMUFiring,2); % length of spike train in samp
text(1,4,plottext)
xlabel('meanDR (pps)')
ylabel('CoVISI (%)')


% List of variables saved
WinLen = LW*fsamp; % window length
%COHT % Coherence function
%F % frequencies
%fsamp % sampling frequency
%Iter % number of iterations
data.MUNumber = size(cutMUFiring,1); % number of motor units
data.MULength = size(cutMUFiring,2); % length of spike train in samp
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
    data.F = m1F;
    data.COF = coherCI;
    data.COF2z = NoiseCI;
    data.Z = m1Z;
    data.cstF = TAcstF;
    data.cstZ = m1CSTz;
    data.deltaband = freqband(1);
    data.alphaband = freqband(2);
    data.betaband = freqband(3);
    data.torqueCorr = torqueCorr;
end

save(filename_save,'data','fsamp', 'WinLen','a','torqueCOV');
print (figSPT(101),'-dpdf',[filename_save(1:end-4) '.pdf']);
end

