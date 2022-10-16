% File created by Matt Topley on 1/29/2022
% This script has been modified to loop through files in directory and process multiple time periods for each file and skip those that are bad trials
% Error and warning handling have also been incorporated in the function
% and nested function. The script also creates pdf directories for each ramp in a common location, een if there are ramp subfolders.
% The actual analysis has also been modified in the following:
% The units are first plotted to better determine contractio time points and
% removal of units, The units recruitment and derecruitment thresholds are set to
% >= 400ms and is attempted at both 5 and 10 pulses, Units firing for less
% than 4s during the time period selected are automatically removed,
% Users can compare cocontraction of antagonists with the
% the actual torque/torque feedback and agonist EMG feedback in the summary figure
% The code has also been modified for simpler variable naming and debugging


function [rere, dere] = Functional_DeltaF_Analysis(fullFileName, usermaxTorque, Auto, saveexcel, keepPics, ramp)
% Subfunction processing delta f values from decomposed MU spike trains.
% Plot MU and delta f comparisons.
% Save results as PDF and Excel files.
%%% v3_CoCo_MT - added comments throughout
%%% Also added antagonist EMG to identify cocontraction
%%% Fixed edge trimming to 400ms and remove trains < 4s long
Pics = keepPics(2);
keepPics = keepPics(1);

close all;

%%% Set up directories - change if > 2 subdirectories
[fileDir,filename,clusterFile, pdfdir] = setupDirectories(fullFileName, 'ramp');

% Open Decomposed MU File
load(fullfile(fileDir,filename),'MUPulses', 'TraceFeedback')
%load(filename, 'MUPulses', 'Torque', 'TorqueFeedback', 'EMGFeedback');
Tableaumap = Tableau;   % colormap if user doesn't have predefined colormap

try
    load(fullfile(fileDir,clusterFile), 'fsamp', 'Torque');
    temp = load(fullfile(fileDir,clusterFile), 'JR3Mz');
    TorqueFeedback = temp.JR3Mz;
end

try
    EMGFeedback = [];
    load(fullfile(fileDir,clusterFile), 'EMG');
end

absRawEMG = abs(EMG(27,:));
rectifiedEMG = absRawEMG-mean(absRawEMG(1:fsamp));
rmsEMG = rms(rectifiedEMG, 500);
EMGFeedback = movmean(rmsEMG, fsamp);
EMGFeedback = EMGFeedback-mean(EMGFeedback(1:fsamp));
fs = 1/fsamp;
[b,a] = butter(3,fs,'low');
EMGthreshold = filtfilt(b,a,rmsEMG);

%%%%% load agonist units
[MUFiring] = SortUnits(MUPulses);

%%%%%
if exist('fsamp') == 0
    fsamp = 2048;
end

checkTorque = Torque-mean(Torque(1:fsamp));
if abs(mean(checkTorque)./usermaxTorque) < 2
    coco = 1;
else
    coco = 0;
end

% Find antagonist file, adjust these if analyzing different muscles
if coco == 1 && contains(filename,'TA')
    [AntagEMGFeedback] = findAntagonistEMG(fullfile(fileDir,clusterFile), 'TA', 'Sol');
        agonistEMGColor = 1;
        antagonistEMGColor = 7;
elseif coco == 1 && contains(filename,'Sol')
     [AntagEMGFeedback] = findAntagonistEMG(fullfile(fileDir,clusterFile), 'Sol', 'TA');
         agonistEMGColor = 7;
         antagonistEMGColor = 1;
     elseif coco == 1 && contains(filename,'MG')
     [AntagEMGFeedback] = findAntagonistEMG(fullfile(fileDir,clusterFile), 'MG', 'TA'); 
         agonistEMGColor = 5;
         antagonistEMGColor = 1;
end

%%%%%
if Auto == 1
[rampStart] = findRamps(TraceFeedback,EMGFeedback,EMG, coco);
if rampStart(1) < 1
    rampStart(1) = 1;
end

%%%% Auto find ramps Comment out if not using auto selection
if ramp == 1
    Trialcut = [rampStart(1) rampStart(1)+30];
elseif ramp == 2
    Trialcut = [rampStart(2) rampStart(2)+30];
else
    Trialcut = [rampStart(3) rampStart(3)+30];
end
stax = Trialcut(1);
disp(Trialcut)
end 

%%%%
% plot units and choose ramp time range in sec

for MU = 1:size(MUFiring,2)
    index = MUFiring{MU};
    figfiring(MU,index) = 1;
end
figfiring(:,length(figfiring):length(EMGFeedback))=0;

f1 = figure(1);
% smallWindow = get(f1, 'Position');
% set(f1, 'Position', get(0,'Screensize'));
for i = 1:length(MUFiring)
    hold all
    plot (MUFiring{i}(2:end)./fsamp,1./diff(MUFiring{i})*fsamp-(20*(i-1)),'.', 'MarkerSize', 10);
    plot([0:length(figfiring)-1]/fsamp,2*fftfilt(hanning(fsamp),figfiring(i,:)')-(20*(i-1)))
    %  text(stax+1, -20*(i-1)+2,num2str(i))
end
hold off
title(clusterFile)

%%% Option to skip current file unless reprocessing another ramp
if Auto ~= 1
%     if ramp == 1
%         skip = input('skip trial? Yes or No: ','s');
%         if contains(skip,'y');
%             return
%         end
%     end
    
    %%% start and end of Ramp
    
    Trialcut = input('Enter start and end of Contraction time (e.g., [5 30]): ');%[60 85];
    if length(Trialcut) ~= 2
        Trialcut = input('Missing time point, Enter start and end of Contraction time (e.g., [5 30]): ');%[60 85];
    end
    
    if Trialcut(1)<10
        stax = 1;
        %     xlim([stax Trialcut(2)+2]);
    else
        stax = Trialcut(1)-2;
        %     xlim([stax Trialcut(2)+2]);
    end
    
    
    for i = 1:length(MUFiring)
        text(stax+1, -20*(i-1)+2,num2str(i))
    end
end

%%%% Remove bad units, either user input or auto
% if you are reassigning numbers for bad MUs, if Auto=1 then dont worry about it
TRUEMU = [1:length(MUFiring)];
rspike = [];
if Auto == 0
    rspike = input('Enter bad motor units (e.g., [1,2,3]): ');
    %badMU = [];
    %badMU = [badMU, userBadMU];
    
    MUFiring(rspike) = [];
    figfiring(rspike,:) = [];
    TRUEMU(rspike) = [];
    
    % replot units with the bad ones removed
    clf(f1);
    for i = 1:length(MUFiring)
        hold all
        plot (MUFiring{i}(2:end)./fsamp,1./diff(MUFiring{i})*fsamp-(20*(i-1)),'.', 'MarkerSize', 20);
        plot([0:length(figfiring)-1]/fsamp,2*fftfilt(hanning(fsamp),figfiring(i,:)')-(20*(i-1)));
        text(stax+1, -20*(i-1)+2,num2str(i))
    end
    title(clusterFile)
    drawnow
    % xlim([stax Trialcut(2)+10]);
    hold off
    set(f1, 'Position', smallWindow);
end

%%%% plot actual torque, torque feedback, and EMG feedback
try
    f2=figure(2);
    subplot(311)
    plot (1/fsamp:1/fsamp:length(Torque)/fsamp,Torque);
    subplot(312)
    plot (1/fsamp:1/fsamp:length(TorqueFeedback)/fsamp,TorqueFeedback);
    subplot(313)
    plot (1/fsamp:1/fsamp:length(EMGFeedback)/fsamp,EMGFeedback)
catch
end

%%%% cut trialcut length is ok
if Trialcut(2)>(length(EMGFeedback)/fsamp)-1
    Trialcut(2) = (floor(length(EMGFeedback)/fsamp))-1;
end

%%%%
%This is where MUFiring, torque, et al should be cut
%based upon the user inputs of ramp start/stop
MUTime = [Trialcut(1) Trialcut(2)];
% TQTime = [(Trialcut(1)*fsamp)+1 (Trialcut(2)*fsamp)];

% if exist(fullfile(fileDir,[filename(1:end-4),'_',num2str(Trialcut(1)),'_',num2str(Trialcut(2)),'_Ramp1.pdf'])) == 2;
%     disp('Ramp 1 exists');
%     return
% elseif exist(fullfile(fileDir,[filename(1:end-4),'_',num2str(Trialcut(1)),'_',num2str(Trialcut(2)),'_Ramp2.pdf'])) == 2;
%     disp('Ramp 2 exists');
%     return
% elseif exist(fullfile(fileDir,[filename(1:end-4),'_',num2str(Trialcut(1)),'_',num2str(Trialcut(2)),'_Ramp3.pdf'])) == 2;
%     disp('Ramp 3 exists');
%     return
% end

try
    %%% cut real torque and Torque Feedback
    [TQTime, timeTorque, cutTorque, cutTorqueFeedback, pfTQfeedback, dfTQfeedback] = TrimTorque(MUTime,Torque, TorqueFeedback, usermaxTorque, coco);
    
    %%% cut EMG feedback
    EMGFeedback = ((EMGFeedback - mean(EMGFeedback(1:fsamp)))./(abs((max(EMGFeedback)-min(EMGFeedback)))))*20;
    %      EMGFeedback = ((rectifiedEMG - mean(rectifiedEMG(1:fsamp)))./(abs((max(rectifiedEMG)-min(rectifiedEMG)))))*20;
    newEMGFeedback = EMGFeedback(TQTime(1):TQTime(2));
    cutTrace = TraceFeedback(TQTime(1):TQTime(2));
    %     EMGthreshold = EMGthreshold(TQTime(1):TQTime(2));
    %     EMGthreshold = EMGthreshold -mean(EMGthreshold (1:fsamp));
    %     newEMGFeedback = rectifiedEMG(TQTime(1):TQTime(2));
    %%% cut antagnoist EMG feedback if applicable
    try
        AntagEMGFeedback = ((AntagEMGFeedback - mean(AntagEMGFeedback(1:fsamp)))./(abs((max(AntagEMGFeedback)-min(AntagEMGFeedback)))))*10;
        newAntagEMGFeedback = AntagEMGFeedback(TQTime(1):TQTime(2));
        
    end
catch
end

maxTorque = max(abs(cutTorqueFeedback));
maxEMG = max(newEMGFeedback);
if coco == 1
    torqueRMSE = sqrt((0-cutTorqueFeedback).^2);
else
    torqueRMSE = sqrt((cutTrace-cutTorqueFeedback).^2);
end
normTorqueRMSE = mean(torqueRMSE)./(max(torqueRMSE)-min(torqueRMSE));

%%% Mean torque of begginning and end of ramps
secondBase = [mean(cutTorque(1:fsamp)),mean(cutTorque((end-fsamp):end))];

%%% Remove spike trains firing for less than 4s (allows for 2s up and 2s down)
[MUFiring, AVERAGE,firing, TRUEMU] = TrimEdges2Ramp(MUTime, TRUEMU, MUFiring);
%%%
% filtered discharge rates and times
rere=[];
dere =[];
for i = 1:length(MUFiring)
    rere(i) = EMGFeedback(MUFiring{i}(1)*fsamp);
    dere(i) = EMGFeedback(MUFiring{i}(end)*fsamp);
    windowDR{i} = AVERAGE(i,min(MUFiring{i}*fsamp):max(MUFiring{i}*fsamp));
    windowST{i} = min(MUFiring{i}):1/fsamp:max(MUFiring{i});
    %     unitvariance(i) = var(windowDR{i});
    %windowTQ{i} = newTorque(MUFiring{i}(1)*fsamp-(TQTime(1)-1):MUFiring{i}(end)*fsamp-(TQTime(1)-1));
end

%%% same figure as the top summary figure plot
if Pics == 4
    
    f3 = figure(3); hold all
    xlim([MUTime(1) MUTime(2)]);
    ylim([0 25]);
    c = 1;
    Tableaumap = Tableau;
    for i = 1:length(windowDR)
        try
            plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:));
            c=c+1;
        catch
            c=1;
            plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:));
        end
        %title({[filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],['All ',num2str(length(MUFiring)),' units']},'interpreter','none')
    end
    %plot(timeTorque,newTorque/max(newTorque)*max(AVERAGE(size(AVERAGE,1),:)),'k','LineWidth',1.75)
    plot(timeTorque,cutTorque,'k','LineWidth',1.75)
    plot([1/fsamp:1/fsamp:size(AVERAGE,2)/fsamp],AVERAGE(size(AVERAGE,1),:),'Color', [.1 .1 .1],'LineWidth',1.75)
    title4All = [filename(1:end-4),'_allunits_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'.pdf'];
    
    if keepPics == 4
        %         saveas(f3,title4All)
        print(f3,'-dpdf',fullfile(pdfdir, title4All))
    end
end

%%%
windowDR{size(AVERAGE,1)} = AVERAGE(size(AVERAGE,1),:);
windowTQ{size(AVERAGE,1)} = secondBase; %%% seems odd
windowST{size(AVERAGE,1)} = 1/fsamp:1/fsamp:(size(AVERAGE,2))/fsamp;

%%%%
% Plot low vs high threshold motor units, find delta f values
try
    MUFiring{size(AVERAGE,1)} = sort(cat(1, MUFiring{:}));
catch
    MUFiring{size(AVERAGE,1)} = sort(cat(2, MUFiring{:}));
end

c=1;
for i = 1:length(windowST)
    for j = 1:length(windowST)
        
        if i == j
            continue
        end
        
        Tdiff = min(windowST{j}) - min(windowST{i});
        if Tdiff>0 % && min(windowST{j})<max(windowST{i})
            % try
            if Pics >= 3
                
                f4 = figure(4);clf; hold all
                try
                    plot(timeTorque,cutTorqueFeedback,'k')
                catch
                    plot(timeTorque,cutTorque,'k')
                end
                plot(windowST{i},windowDR{i},'b','linewidth',2);hold all
                try
                    plot(MUFiring{i}(2:end),1./(diff(MUFiring{i})),'b.')
                end
                plot(windowST{j},windowDR{j}+10,'r','linewidth',2)
                try
                    plot(MUFiring{j}(2:end),1./(diff(MUFiring{j}))+10,'r.')
                end
                xlim([MUTime(1) MUTime(2)])
                ylim([0 30])
            end
            
            Fstart = windowDR{i}(round(Tdiff*fsamp));
            Tstart = windowST{i}(windowDR{i}==Fstart);
            
            Fmax = max(windowDR{i}(round(Tdiff*fsamp):end));
            Tmax = windowST{i}(windowDR{i}==Fmax);
            
            ratemod = Fmax - Fstart;
            
            if max(windowST{j})<max(windowST{i})
                Tend = max(windowST{j});
            else
                Tend = max(windowST{i});
            end
            
            Fend = windowDR{i}(round(windowST{i}*fsamp)==round(Tend*fsamp));
            deltaf =  Fstart - Fend;
            
            plot([Tstart windowST{j}(1)],[Fstart, windowDR{j}(1)+10],'k--')
            plot([Tend windowST{j}(end)],[Fend, windowDR{j}(end)+10],'k--')
            plot([Tstart windowST{i}(end)+1],[Fstart Fstart],'k:')
            plot([Tstart windowST{i}(end)+1],[Fend Fend],'k:')
            plot([Tmax Tmax],[Fmax Fstart],'k','linewidth',2)
            
            ratemodtxt = sprintf('%0.1f pps',ratemod);
            text(Tmax+1,Fmax+1,ratemodtxt)
            plot([Tend+1 windowST{j}(end)+1],[Fend, Fstart],'k','linewidth',2)
            deltaftxt = sprintf('%0.1f pps',deltaf);
            text(Tend+1.2,(Fend+Fstart)/2,deltaftxt)
            plot([windowST{i}(1) windowST{j}(1)],[windowDR{j}(1)+10 windowDR{j}(1)+10],'k','linewidth',2)
            Tdifftxt = sprintf('%0.1f s',Tdiff);
            text(windowST{i}(1),windowDR{j}(1)+9,Tdifftxt)
            
            initiali = windowDR{i}(1);
            maxi = max(windowDR{i});
            finali = windowDR{i}(end);
            initialj = windowDR{j}(1);
            maxj = max(windowDR{j});
            finalj = windowDR{j}(end);
            
            try
                
                %%% find peak on CST
                [M,I] = max(AVERAGE(end,:));
                
                %%% Segment windowDRicut
                thresh = I/fsamp;
                
                windowDRicut = windowDR{i}((windowST{i}>=thresh));
                windowDRjcut = windowDR{j}((windowST{j}>=thresh));
                
                windowSTicut = windowST{i}((windowST{i}>=thresh));
                windowSTjcut = windowST{j}((windowST{j}>=thresh));
                
                if windowSTjcut(1)<thresh
                    windowDRicut = windowDRicut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    windowDRjcut = windowDRjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                    windowSTicut = windowSTicut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    windowSTjcut = windowSTjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                else
                    windowDRicut = windowDRicut(round(1+((windowSTjcut(1)-thresh)*fsamp):min([length(windowDRicut) length(windowDRjcut)])+((windowSTjcut(1)-thresh)*fsamp)));
                    windowDRjcut = windowDRjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                    windowSTicut = windowSTicut(round(1+((windowSTjcut(1)-thresh)*fsamp):min([length(windowDRicut) length(windowDRjcut)])+((windowSTjcut(1)-thresh)*fsamp)));
                    windowSTjcut = windowSTjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                end
                
                %%% also part of same figure
                plot(windowSTicut,windowDRicut,'b','linewidth',4)
                plot(windowSTjcut,windowDRjcut+10,'r','linewidth',4)
                %windowDRicut = windowDR{i}(round(Tdiff*fsamp):round((Tend*fsamp)-(min(windowST{i})*fsamp)));
                [rvalue slope inter] = regression(windowDRicut,windowDRjcut);
            catch
                rvalue = 0;
%                 slope = 0;
            end
            
            titletxt1 = sprintf('Units %0.0f and %0.0f - deltaf = %0.1f, rate-rate slope = %0.2f, rate-rate r-value = %0.2f',TRUEMU(i),TRUEMU(j),deltaf,slope,rvalue);
            titletxt2 = sprintf('Unit %0.0f - Tstart = %0.1f, Tend = %0.1f, Fstart = %0.1f, Fmax = %0.1f, Fend = %0.1f',TRUEMU(i),windowST{i}(1),windowST{i}(end),initiali,maxi,finali);
            titletxt3 = sprintf('Unit %0.0f - Tstart = %0.1f, Tend = %0.1f, Fstart = %0.1f, Fmax = %0.1f, Fend = %0.1f',TRUEMU(j),windowST{j}(1),windowST{j}(end),initialj,maxj,finalj);
            title({titletxt1,titletxt2,titletxt3})
            
            
            %%% Output of variables for excel file
            if or(Tdiff > 1-2 && rvalue > 0.836-2 && ratemod > 0.5-2 && i~=length(windowST),i==length(windowST)); %%%%%% changed
                %                 output(c,:) = [TRUEMU(i),TRUEMU(j),deltaf,Tdiff,rvalue,slope,ratemod,initiali,maxi,finali,initialj,maxj,finalj,Fstart,...
                %                     min(windowST{i}),max(windowST{i}),min(windowST{j}),max(windowST{j}),...
                %                     newTorque((min(MUFiring{i})*fsamp-TQTime(1))),newTorque((max(MUFiring{i})*fsamp-TQTime(1))),...
                %                     newTorque((min(MUFiring{j})*fsamp-TQTime(1))),newTorque((max(MUFiring{j})*fsamp-TQTime(1))),...
                %                     newEMGFeedback((min(MUFiring{i})*fsamp-TQTime(1))),newEMGFeedback((max(MUFiring{i})*fsamp-TQTime(1))),...
                %                     newEMGFeedback(min(MUFiring{j})*fsamp-TQTime(1)),newEMGFeedback((max(MUFiring{j})*fsamp-TQTime(1)))];%%%%%%%%%%
                output(c,:) = [TRUEMU(i),TRUEMU(j),deltaf,Tdiff,rvalue,slope,ratemod,initiali,maxi,finali,initialj,maxj,finalj,Fstart,...
                    min(windowST{i}),max(windowST{i}),min(windowST{j}),max(windowST{j}),...
                    cutTorque((min(MUFiring{i})*fsamp-TQTime(1))),cutTorque((max(MUFiring{i})*fsamp-TQTime(1))),...
                    cutTorque((min(MUFiring{j})*fsamp-TQTime(1))),cutTorque((max(MUFiring{j})*fsamp-TQTime(1))),...
                    newEMGFeedback((min(MUFiring{i})*fsamp-TQTime(1))),newEMGFeedback((max(MUFiring{i})*fsamp-TQTime(1))),...
                    newEMGFeedback(min(MUFiring{j})*fsamp-TQTime(1)),newEMGFeedback((max(MUFiring{j})*fsamp-TQTime(1)))];%%%%%%%%%%
                c=c+1;
            end
            
            if keepPics >= 3
                savename = [filename(1:end-4),'_Unit_', num2str(TRUEMU(i), '%02i'),'_Unit_',num2str(TRUEMU(j), '%02i'),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))];
                print(f4,'-dpdf',fullfile(pdfdir, savename))
                
                %                 print([filename(1:end-4),'_Unit_', num2str(TRUEMU(i), '%02i'),'_Unit_',num2str(TRUEMU(j), '%02i'),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],'-dpdf')
                %
                %                 saveas (f4,[savename,'.png'])
            end
        end
    end
end

if Pics >= 2 % IF u want PDFs of individual units displayed
    
    plot_axis = [MUTime(1) MUTime(2) 0 20];
    axis(plot_axis)
    %savename = [filename(1:end-4),'_all units_win'];
    c=1;
    for i = 1:length(MUFiring)
        starttime = MUFiring{i}(1);
        endtime = MUFiring{i}(end);
        f99 = figure(99);clf;
        
        if c<21
            plot(MUFiring{i}(2:end),1./(diff(MUFiring{i})),'o','Color',Tableaumap(c,:))
            hold on;
            plot(starttime:1/fsamp:endtime,AVERAGE(i,starttime*fsamp:endtime*fsamp),'Color',Tableaumap(c,:),'LineWidth', 2)
        else
            c=1;
            plot(MUFiring{i}(2:end),1./(diff(MUFiring{i})),'o','Color',Tableaumap(c,:))
            hold on;
            plot(starttime:1/fsamp:endtime,AVERAGE(i,starttime*fsamp:endtime*fsamp),'Color',Tableaumap(c,:),'LineWidth', 2)
        end
        title({filename(1:end-4), ['Unit ', num2str(TRUEMU(i), '%02i')]}, 'Interpreter', 'none');
        axis(plot_axis)
        
        if keepPics >= 2
            savenamef99 = [filename(1:end-4),'_Unit_', num2str(TRUEMU(i), '%02i'),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))];
            print(f99,'-dpdf',fullfile(pdfdir, savenamef99))
        end
        c=c+1;
    end
end

%%% uncertain of purpose of this section as of 1/29/22
[M,I] = max(cutTorqueFeedback(5*fsamp:end));
I= I+5*fsamp;
% plot(timeTorque(I-9*fsamp:I-fsamp),detrend(newTorque(I-9*fsamp:I-fsamp))); hold all
% plot(timeTorque(I-9*fsamp:I-fsamp),detrend(newTorque(I+fsamp:I+9*fsamp)))
try
    ascendingSD = std(detrend(cutTorque(I-10*fsamp:I)));
    descendingSD = std(detrend(cutTorque(I:I+10*fsamp)));
catch
    ascendingSD = 0;
    descendingSD = 0;
end
try
    A = cutTorque(I-10*fsamp:I);
    D = cutTorque(I:I+10*fsamp);
    [ascendingFIT,m1,b1] = regression(1:length(A)',A);
    [decendingFIT,m2,b2] = regression(1:length(D)',D);
catch
    ascendingFIT = 0;
    decendingFIT = 0;
end
[maxTorque,IT] = max(cutTorque);

peakdiff = (IT-I)/fsamp;

Torqueout = zeros(1,size(output,2));
Torqueout(1) = 0;
Torqueout(2) = 0;
Torqueout(3) = maxTorque;
Torqueout(4) = peakdiff;
Torqueout(5) = ascendingSD;
Torqueout(6) = descendingSD;
Torqueout(7) = ascendingFIT;
Torqueout(8) = decendingFIT;
Torqueout(9) = usermaxTorque;
Torqueout(10) = max(newEMGFeedback);

output(size(output,1)+1,:) = Torqueout;

%%% plot and save sumary figure
if Pics >= 1
    sumFig = figure(101); hold all
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1) pos(2)-pos(4) pos(3) pos(4)*2])
    subplot(8,2,[1 4]); hold all
    xlim([MUTime(1) MUTime(2)]);
    if coco == 1
        ylim([-5 25]);
    else
        ylim([0 25]);
    end
    c = 1;
    Tableaumap = Tableau;
    
    for i = 1:length(windowDR)
        try
            plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:));
            c=c+1;
        catch
            c=1;
            plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:));
        end
        title({[filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],['All ',num2str(length(MUFiring)),' units']},'interpreter','none')
    end
    
    try
        %         plot(timeTorque,newTorqueFeedback,'MarkerColor', TQlinecolor,'LineWidth',2)%
        plot(timeTorque,dfTQfeedback,'k','LineWidth',2)%
        plot(timeTorque,pfTQfeedback,'r','LineWidth',2)%
        plot(timeTorque,newEMGFeedback,'SeriesIndex', agonistEMGColor,'LineWidth',1,'LineStyle',':')
        plot(timeTorque,newAntagEMGFeedback,'SeriesIndex',antagonistEMGColor,'LineStyle',':','LineWidth',1)
    catch
        %         plot(timeTorque,newTorqueFeedback,'Color',TQlinecolor,'LineWidth',2)
        plot(timeTorque,dfTQfeedback,'k','LineWidth',2)%
        plot(timeTorque,pfTQfeedback,'r','LineWidth',2)%
        plot(timeTorque,newEMGFeedback,'Color', 'k','LineWidth',1.5,'LineStyle',':')
    end
    plot([1/fsamp:1/fsamp:size(AVERAGE,2)/fsamp],AVERAGE(size(AVERAGE,1),:),'Color', [.1 .1 .1],'LineWidth',1.75)
    
    title4All = [filename(1:end-4),'_allunits_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'.pdf'];
    title(title4All)
    
    alloutput = output;
    outputclean = alloutput;
    outputclean(outputclean(:,1) == 0,:) = [];
%     outputAverage = outputclean;
%     outputAverage(outputAverage(:,4)<1,:) = [];
%     outputAverage(outputAverage(:,16)-outputAverage(:,18)<1.5,:) = [];
    
    %%%%
    % Replacing all unit comparisons with the average delta f for each unit
    avgdfs = [];
    Units = unique(outputclean(:,2));
    for jk = 1:length(Units)
        placeholder = [];
        placeholder = outputclean(outputclean(:,2) == Units(jk),:);
        if size(placeholder,1) < 2
            avgdfs(jk,:) = placeholder;
        else
            avgdfs(jk,:) = mean(placeholder);
        end
    end
    
    %%% Plotting delta f histograms
    deltaf_binwidth = [-20:1:20];
    subplot(825); hold all
    bar(deltaf_binwidth,histc(outputclean(:,3),deltaf_binwidth)); %outputclean - unit avgs via avgdfs
    figtitle = sprintf(['delta-f  %0.2f ' char(177) ' %0.2f'], [mean(avgdfs(:,3));std(avgdfs(:,3))]);
    title(figtitle)
    xlim([min(deltaf_binwidth) max(deltaf_binwidth)])
    
    deltaf_binwidth = [0:0.2:10];
    subplot(826)
    bar(deltaf_binwidth,histc(outputclean(:,6),deltaf_binwidth));
    figtitle = sprintf(['Rate-rate slope  %0.2f ' char(177) ' %0.2f'], [mean(avgdfs(:,6));std(avgdfs(:,6))]);
    title(figtitle)
    xlim([min(deltaf_binwidth) max(deltaf_binwidth)])
    
    deltaf_binwidth = [0:0.2:10];
    subplot(827)
    bar(deltaf_binwidth,histc(outputclean(:,4),deltaf_binwidth));
    figtitle = sprintf(['Time difference  %0.2f ' char(177) ' %0.2f'], [mean(avgdfs(:,4));std(avgdfs(:,4))]);
    title(figtitle)
    xlim([min(deltaf_binwidth) max(deltaf_binwidth)])
    
    %%% Discharge rates output unmodified
    data_MU = alloutput;
    data_MU(data_MU(:,1) ~= 0,:) = [];
    data_MU(data_MU(:,2) == 0,:) = [];
    
    deltaf_binwidth = [0:0.5:30];
    subplot(828)
    bar(deltaf_binwidth,histc(data_MU(:,12),deltaf_binwidth))
    figtitle = sprintf(['Max discharge rate  %0.2f ' char(177) ' %0.2f'], [mean(data_MU(:,12));std(data_MU(:,12))]);
    title(figtitle)
    xlim([min(deltaf_binwidth) max(deltaf_binwidth)])
    
    %%% Adjusting the figures sizing
    subplot(8,2,7)
    subpos = get(gca,'position');          % gca points at the second one
    subpos(2) = subpos(2)-subpos(4)*6;              % reduce the height by half
    subpos(4) = subpos(4)*2;              % reduce the height by half
    set(gca,'position',subpos);
    
    subplot(8,2,8)
    subpos = get(gca,'position');          % gca points at the second one
    subpos(2) = subpos(2)-subpos(4)*6;              % reduce the height by half
    subpos(4) = subpos(4)*2;              % reduce the height by half
    set(gca,'position',subpos);
    
    subplot(8,2,5)
    subpos = get(gca,'position');          % gca points at the second one
    subpos(2) = subpos(2)-subpos(4)*4;              % reduce the height by half
    subpos(4) = subpos(4)*2;              % reduce the height by half
    set(gca,'position',subpos);
    
    subplot(8,2,6)
    subpos = get(gca,'position');          % gca points at the second one
    subpos(2) = subpos(2)-subpos(4)*4;              % reduce the height by half
    subpos(4) = subpos(4)*2;              % reduce the height by half
    set(gca,'position',subpos);
    
    subplot(8,2,[1 4])
    subpos = get(gca,'position');          % gca points at the second one
    subpos(2) = subpos(2)-subpos(4);              % reduce the height by half
    subpos(4) = subpos(4)*2;              % reduce the height by half
    set(gca,'position',subpos);
    
    if keepPics >= 1
        
        if ramp == 1 
            sumFigTitle = [filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp1.pdf'];
        elseif ramp == 2
            sumFigTitle = [filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp2.pdf'];
        elseif ramp == 3
            sumFigTitle = [filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp3.pdf'];
        end
        print(sumFig,'-dpdf',fullfile(fileDir, sumFigTitle))
        %saveas(sumFig,title4All)
    end
end

%%% Save the excel file
if saveexcel == 1
        exceldir = [fileDir '\results_excel'];
    if ~exist([fileDir  '\results_excel'],'dir')
        mkdir([fileDir '\results_excel']);
    end

    
    %     if MUTime(1) < 20
    dlmwrite(fullfile(exceldir,[filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_deltaf.xls']),output,'delimiter','\t')
    %         elseif MUTime(1) < 60
    %             sumFigTitle = [filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp2.pdf'];
    %         else
    %             sumFigTitle = [filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp3.pdf'];
    %         end
    
    %     try
    %         exdir = [Feedbackdir '/results_excel'];
    %         try
    %             cd (exdir)
    %             dlmwrite([filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_deltaf.xls'],output,'delimiter','\t')
    %             cd (currentdir)
    %         catch
    %             mkdir('results_excel')
    %             cd (exdir)
    %             dlmwrite([filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_deltaf.xls'],output,'delimiter','\t')
    %             cd (currentdir)
    %         end
    %     end
end

% Load colormap if you don't have start-up file
    function [tableau] = Tableau
        tableau =     [0.1216    0.4667    0.7059
            0.6824    0.7804    0.9098
            1.0000    0.4980    0.0549
            1.0000    0.7333    0.4706
            0.1725    0.6275    0.1725
            0.5961    0.8745    0.5412
            0.8392    0.1529    0.1569
            1.0000    0.5961    0.5882
            0.5804    0.4039    0.7412
            0.7725    0.6902    0.8353
            0.5490    0.3373    0.2941
            0.7686    0.6118    0.5804
            0.8902    0.4667    0.7608
            0.9686    0.7137    0.8235
            0.4980    0.4980    0.4980
            0.7804    0.7804    0.7804
            0.7373    0.7412    0.1333
            0.8588    0.8588    0.5529
            0.0902    0.7451    0.8118
            0.6196    0.8549    0.8980];
    end
end

