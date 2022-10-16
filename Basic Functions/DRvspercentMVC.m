function [] = Batch_Limbslopes
clear all;% -except c %%close all
filenames = ls('*.xls')
%filenames = ls('CoCoTest09_44426_RTA_1_v23decomposed_MUCLEANED.mat')
Tableaumap = Tableau;
allLimbs = [];
%c = input('Number? ');
for i = 1:size(filenames,1)
    try
        [slops] = DRvspercentMVC(filenames(i,:))
        allLimbs = [allLimbs;slops];
    end 
end
    fig_binwidth = [-20:1:20];
    hold all
    legend('-DynamicLegend', 'Location', 'bestoutside');
    % legName = (['All Dots Ramps, Total Units: ', num2str(totUnit)]) %filenames(i,1:11);
    % nums = sprintf([' delta-f  %0.2f ' char(177) ' %0.2f'], [mean(allDFs(:,3));std(allDFs(:,3))]);
    % dname = [legName nums];
    bar(fig_binwidth,histc(abs(allLimbs(:,1)),fig_binwidth), 'FaceAlpha', 0.5, 'DisplayName', 'Limb1' ); 
    bar(fig_binwidth,histc(abs(allLimbs(:,2)),fig_binwidth), 'FaceAlpha', 0.5, 'DisplayName', 'Limb2');
    bar(fig_binwidth,histc(abs(allLimbs(:,3)),fig_binwidth), 'FaceAlpha', 0.5, 'DisplayName', 'Limb3');
    bar(fig_binwidth,histc(abs(allLimbs(:,4)),fig_binwidth), 'FaceAlpha', 0.5, 'DisplayName', 'Limb4');
    %title('All Subject CoCo vs EMG Ramps')
    xlabel('Average Pairwise Delta Fs')
    ylabel('Count')
    xlim([min(fig_binwidth) max(fig_binwidth)])
    hold off
end

function [slops] = DRvspercentMVC(filename)
% window = 0;
output = load(filename);
currdir = cd;
decompDir = fileparts(currdir);

Under = strfind(filename,'_');
% DecompFile = [filename(1:(Under(6)-1)),'.mat'];
startRamp = str2num(filename(Under(6)+1:Under(7)-1));
endRamp = str2num(filename(Under(7)+1:Under(8)-1));

load(fullfile(decompDir, DecompFile));
Trialcut = [startRamp endRamp];

% clear all; close all;
% filename = uigetfile;
load(filename)

Under = strfind(filename,'_');

clusFile = [filename(1:(Under(4)-1)),'.mat'];
load(clusFile);

% Trialcut = input('Contraction time: ')
% Trialcut = [Trialcut(1), Trialcut(2)];
% MUTime = [Trialcut(1) Trialcut(2)];
% TQTime = [(Trialcut(1)*fsamp) (Trialcut(2)+1)*fsamp];
MUTime = [Trialcut(1) Trialcut(2)];
TQTime = [(Trialcut(1)*fsamp)+1 Trialcut(2)*fsamp];

cutTorque = JR3Mz(TQTime(1):TQTime(2));
newTorque = (cutTorque-(mean(cutTorque(1:fsamp))));
timeTorque = TQTime(1)/fsamp:1/fsamp:MUTime(2);
if mean(newTorque)<0
    newTorque = -newTorque;
end
% mvc = -(JR3Mz-(mean(JR3Mz(1:fsamp))));
% mvc= (mvc./max(mvc));

goodEMG = std(EMG, 0, 2);
goodEMG = EMG(goodEMG<100,:);

%Filtering agonist EMG
sumdifEMG = sum(abs(goodEMG(:,[min(TQTime):max(TQTime)])));
sumdifEMG = sumdifEMG-mean(sumdifEMG(1:fsamp));
sumdifEMG = sumdifEMG/max(sumdifEMG)*10;
[b,a]   = butter(3,1/fsamp/2);
sumEMG = filtfilt(b,a,sumdifEMG);
sumEMG = sumEMG/max(sumEMG)*10;
% 
% figure;
% tiledlayout(2,1)
% nexttile(1)
% plot(timeTorque, newTorque)
% nexttile(2)
% plot(timeTorque, sumEMG)
% hold on
% plot(timeTorque,sumdifEMG);
% hold off

fs = 1/fsamp;
[~,DRlen] = sort(cellfun(@length,MUPulses));
MUFiring = MUPulses(DRlen);
TRUEMU = [1:length(MUFiring)];

%timeTorque = TQTime(1)/fsamp:1/fsamp:MUTime(2);
% ISI = cellfun(@(x) diff(x/fsamp), MUPulses2, 'UniformOutput', false);
% stdISI = cellfun(@std, ISI, 'UniformOutput', false);
% meanISI = cellfun(@mean, ISI, 'UniformOutput', false);
% IDR = cellfun(@(x) 1./(diff(x)/fsamp), MUPulses2, 'UniformOutput', false);
% stdIDR = cellfun(@std, IDR, 'UniformOutput', false);
% meanIDR = cellfun(@mean, IDR, 'UniformOutput', false);
% timeMUPulses = cellfun(@(x) x/fsamp, MUPulses2, 'UniformOutput', false);
for i = 1:length(MUFiring)
    loopMuF{i} = MUFiring{i}(MUFiring{i}>Trialcut(1).*fsamp & MUFiring{i}<Trialcut(2).*fsamp);
end

removetrains = cellfun(@length, loopMuF);
% removetrains<10

% TRUEMU = [TRUEMU(~cellfun(@isempty, loopMuF)),0];
% loopMuF = loopMuF(~cellfun(@isempty, loopMuF));

TRUEMU(removetrains<10) = [];
loopMuF(removetrains<10) = [];

for i = 1:length(loopMuF) % Create cell array of MU datapoints in seconds
    
    MuF = loopMuF{i};
    preTSAISI = (diff(MuF/fsamp));
    clear r rr; [rr r] = find(preTSAISI(1:5)>=0.8);
    
    if isempty(rr);
        preMUPulses = MuF;
    else
        y = max(r);
        preMUPulses = MuF(y+1:end);
    end
    
    clear h hh
    [hh h] = find(preTSAISI(end-5:end)>=0.8);
    if isempty(hh);
        preMUPulses = preMUPulses;
    else
        z = min(h);
        preMUPulses = preMUPulses(1:end-5+z-2);
    end
    MuF = preMUPulses;
    loopMUFiring{i} = MuF/fsamp;
    
end

TRUEMU = [TRUEMU(~cellfun(@isempty, loopMUFiring)),0];
MUFiring = loopMUFiring(~cellfun(@isempty, loopMUFiring));

for i = 1:length(MUFiring)
    MuF = MUFiring{i};
    index = round(MuF*fsamp);
    firing(i,index) = 1;
end

firing(:,length(firing)+fsamp)=0; %adds a second of 0s to the end
AVERAGE = 1./sum(hanning(fsamp*2))*filtfilt(hanning(fsamp*2),1,firing')';%AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);

for i = 1:length(MUFiring)
    windowDR{i} = AVERAGE(i,min(MUFiring{i}*fsamp):max(MUFiring{i}*fsamp));
    windowST{i} = min(MUFiring{i}):1/fsamp:max(MUFiring{i});
    %windowTQ{i} = newTorque(MUFiring{i}(1)*fsamp-(TQTime(1)-1):MUFiring{i}(end)*fsamp-(TQTime(1)-1));
end

currdir = pwd;
slops = [];
for jk = 1:length(MUFiring)
    %force = mvc;
    
    try
        [ymax,updown] = max(windowDR{jk});
        half = windowST{jk}(updown);
        
        if half-MUFiring{jk}(2)> 4 && MUFiring{jk}(end)-half > 4
            L1 = MUFiring{jk}(2)+1;
            limb1 = MUFiring{jk}(MUFiring{jk}<L1);
            fitL1 = polyfit(limb1(2:end),1./diff(limb1),1);
            L1Fitline = polyval(fitL1,limb1(2:end));
            
            %L2 = MUFiring{jk}(half)+1;
            limb2 = MUFiring{jk}(MUFiring{jk}<half);
            limb2 = limb2(limb2>=L1);
            fitL2 = polyfit(limb2(2:end),1./diff(limb2),1);
            L2Fitline = polyval(fitL2,limb2(2:end));
            
            limb3 = MUFiring{jk}(MUFiring{jk}>= half);
            L3 = limb3(end)-1;
            limb4 =  limb3(limb3> L3);
            [test, L4] = max(1./diff(limb4));
            L4 = limb4(L4);
            limb4 = limb4(limb4>=L4);
            limb3 =  limb3(limb3< L4);
            
            fitL3 = polyfit(limb3(2:end),1./diff(limb3),1);
            L3Fitline = polyval(fitL3,limb3(2:end));
            fitL4 = polyfit(limb4(2:end),1./diff(limb4),1);
            L4Fitline = polyval(fitL4,limb4(2:end));
            
            p = figure ();
            tiledlayout(1,2)
            
            nexttile(1)
            hold on
            plot(limb1(2:end),1./(diff(limb1)),'o', 'Color', 'k');
            plot(limb1(2:end),L1Fitline, 'Color',  [.5 .5 .5], 'LineWidth', 2);
            plot(limb2(2:end),1./(diff(limb2)),'o', 'Color', 'k');
            plot(limb2(2:end),L2Fitline, 'Color',  [.3 .3 .3], 'LineWidth', 2);
            plot(limb3(2:end),1./diff(limb3),'o', 'Color', 'r');
            plot(limb3(2:end),L3Fitline, 'Color', 'r', 'LineWidth', 2);
            plot(limb4(2:end),1./diff(limb4),'o', 'Color', 'r');
            plot(limb4(2:end),L4Fitline, 'Color', 'r', 'LineWidth', 3);
            hold off
            
            xlim([L1-3 limb4(end)+2]);
            ylim([0 ymax+5]);
            ylabel('Instantaneous Discharge Rate');
            xlabel('Time');
            
            nexttile(2)
            hold on
            %         plot(sumEMG.*fsamp(MUFiring{jk}(2:half)./fsamp),1./(diff(MUFiring{jk}(1:half))),'o', 'Color', 'k');
            %         plot(sumEMG(MUFiring{jk}((half)+1:end).*fsamp)./1000,1./(diff(MUFiring{jk}(half:end))),'o', 'Color', 'r');
            plot(sumEMG(limb1(2:end)*fsamp),1./diff(limb1), 'o', 'Color', 'k')
            plot(sumEMG(limb2(2:end)*fsamp),1./diff(limb2(1:end)), 'o', 'Color',  'k')
            plot(sumEMG(limb3(2:end)*fsamp),1./diff(limb3(1:end)), 'o', 'Color',  'r')
            plot(sumEMG(limb4(2:end)*fsamp),1./diff(limb4(1:end)), 'o', 'Color',  'r')
            %         plot(sumEMG(MUFiring{jk}(half+1:end)*fsamp),1./diff(MUFiring{jk}(half:end)), 'o', 'r')
            %         plot(sumEMG(1:half),1./(diff(MUFiring{jk}(1:half))),'o', 'Color', 'k');
            %         plot(sumEMG(half:end),1./(diff(MUFiring{jk}(half:end))),'o', 'Color', 'r');
            hold off
            xcut = max(sumEMG);
            xlim([-1 xcut+1]);
            ylim([0 ymax+5]);
            ylabel('Instantaneous Discharge Rate');
            xlabel('EMG Trace');
           
            %         nexttile(3)
            %         hold on
            %         plot(mvc(MUFiring{jk}(2:half)*fsamp).*10,1./(diff(MUFiring{jk}(1:half))),'o', 'Color', 'k');
            %         plot(mvc(MUFiring{jk}((half)+1:end)*fsamp).*10,1./(diff(MUFiring{jk}(half:end))),'o', 'Color', 'r');
            %         hold off
            %         xstart = min(mvc);
            %         xend = max(mvc)*10;
            %
            %         xlim([xstart-1 xend]);
            %         ylim([0 ymax]);
            %         ylabel('Instantaneous Discharge Rate');
            %         xlabel('Percentage of MVC');
            %title4All = ['UpDownLimbs_Unit_',num2str(jk),'_','.pdf']
        end
         slops(jk,:)=[fitL1(1),fitL2(1),fitL3(1),fitL4(1)];
    catch
        continue
    end
    slops(slops(:,1) == 0,:) = [];
    %      cd PDFs
    %      saveas(p,title4All)
    %     cd (currdir)
end

end


