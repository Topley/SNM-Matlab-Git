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


function [deltaFData] = delta_f_analysis_TRD_v1(fullFilepath, usermaxTorque, Auto, saveexcel, keepPics)
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
[kkFile, clusterFile, decompFile, cleanFile, ~] = Get_TrialFiles(fullFilepath, 'ramp');
[fileDir, fullFileName] = fileparts(decompFile);
under = strfind(fullFileName, '_');
SubjectID = fullFileName(1:under(2)-1);
% Open Decomposed MU File
cleanVars = matfile(cleanFile);
MUPulses = cleanVars.MUPulses;

%%%%%
if exist('fsamp') == 0
    fsamp = 2048;
end

try
    clusterVars = matfile(clusterFile);
    EMG = clusterVars.EMG;
    [~, EMGFeedback] = EMGpro(EMG, 'filter', {'butter', 4, 5/1024, 'low'});
end

kkVars = matfile(kkFile);
copS = kkVars.COP;
cpAP = copS.WeightedY;
unpadAP = cpAP ~= 0;
paddAP = cpAP == 0;
offset = min(cpAP(unpadAP));
cpAP(paddAP) = offset;
cpAP = cpAP - mean(cpAP(1:fsamp));
TorqueFeedback = cpAP;

%%%%% load agonist units
[MUFiring] = SortUnits(MUPulses);

%%%%%
if Auto == 1
    [Trialcut] = Set_RampWindow(TorqueFeedback, EMGFeedback, EMG, 0);
    
    if Trialcut(1) < 1
        Trialcut(1) = 1;
    end
    totalRamps = length(Trialcut);
else
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
    endax = Trialcut(2);
end

if exist('totalRamps','var' )
    %%%% Auto find ramps Comment out if not using auto selection
    for jj = 1:totalRamps
        
        try
            rampStart(jj) = Trialcut(jj);
            rampEnd(jj) = rampStart(jj+1) - 5;
        catch
            rampEnd(jj) = Trialcut(jj) + 30;
        end
    end
end

%%%% Remove bad units, either user input or auto
% if you are reassigning numbers for bad MUs, if Auto=1 then dont worry about it

if Auto == 0
    TRUEMU = [1:length(MUFiring)];
    rspike = [];
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
    title(SubjectID)
    drawnow
    % xlim([stax Trialcut(2)+10]);
    hold off
    %set(f1, 'Position', smallWindow);
end

% %%%% cut trialcut length is ok
% if Trialcut(2)>(length(EMGFeedback)/fsamp)-1
%     Trialcut(2) = (floor(length(EMGFeedback)/fsamp))-1;
% end
Cmap = Tableau;
%This is where MUFiring, torque, et al should be cut
%based upon the user inputs of ramp start/stop
for jj = 1 : totalRamps
    stax = rampStart(jj);
    endax = rampEnd(jj);
    MUTime = [stax endax];
    
    TRUEMU = [1:length(MUFiring)];
    rspike = [];
    
    auxTimeCut = MUTime .* fsamp;
    
    %%% cut real torque and Torque Feedback
    [cutTorque,timeTorque] = TrimAuxChannel(MUTime, TorqueFeedback);
    [cutEMG,~] = TrimAuxChannel(MUTime, EMGFeedback);
    % EMGFeedback = ((EMGFeedback - mean(EMGFeedback(1:fsamp)))./(abs((max(EMGFeedback)-min(EMGFeedback)))))*20;
    % newEMGFeedback = EMGFeedback(auxTimeCut(1):auxTimeCut(2));
    
    %%% Mean torque of begginning and end of ramps
    secondBase = [mean(cutTorque(1:fsamp)),mean(cutTorque((end-fsamp):end))];
    
    %%% Remove spike trains firing for less than 4s (allows for 2s up and 2s down)
    [MUFiring2, MUTrain, MUSpikes, TRUEMU] = TrimEdges2Ramp(MUTime, TRUEMU, MUFiring);
    %%%
    windowDR = {};
    windowST = {};
    for i = 1:length(MUFiring2)
        
        windowDR{i} = MUTrain(i,min(MUFiring2{i}*fsamp):max(MUFiring2{i}*fsamp));
        windowST{i} = min(MUFiring2{i}):1/fsamp:max(MUFiring2{i});
        
    end
    
    %%% same figure as the top summary figure plot
    if Pics == 4
        
        f3 = figure(3); hold all
        xlim([MUTime(1) MUTime(2)]);
        ylim([0 25]);
        c = 1;
        Cmap = Tableau;
        for i = 1:length(windowDR)
            try
                plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Cmap(c,:));
                c=c+1;
            catch
                c=1;
                plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Cmap(c,:));
            end
        end
        plot(timeTorque,cutTorque,'k','LineWidth',1.75)
        plot([1/fsamp:1/fsamp:size(MUTrain,2)/fsamp],MUTrain(size(MUTrain,1),:),'Color', [.1 .1 .1],'LineWidth',1.75)
        
        if keepPics == 4
            %         saveas(f3,title4All)
            print(f3,'-dpdf',fullfile(pdfdir, title4All))
        end
    end
    
    %%%
    windowDR{size(MUTrain,1)} = MUTrain(size(MUTrain,1),:);
    windowTQ{size(MUTrain,1)} = secondBase; %%% seems odd
    windowST{size(MUTrain,1)} = 1/fsamp:1/fsamp:(size(MUTrain,2))/fsamp;
    
    %%%%
    % Plot low vs high threshold motor units, find delta f values
    try
        MUFiring2{size(MUTrain,1)} = sort(cat(1, MUFiring2{:}));
    catch
        MUFiring2{size(MUTrain,1)} = sort(cat(2, MUFiring2{:}));
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
                        plot(timeTorque,cutTorque,'k')
                    catch
                        plot(timeTorque,cutTorque,'k')
                    end
                    plot(windowST{i},windowDR{i},'b','linewidth',2);hold all
                    try
                        plot(MUFiring2{i}(2:end),1./(diff(MUFiring2{i})),'b.')
                    end
                    plot(windowST{j},windowDR{j}+10,'r','linewidth',2)
                    try
                        plot(MUFiring2{j}(2:end),1./(diff(MUFiring2{j}))+10,'r.')
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
                    [M,I] = max(MUTrain(end,:));
                    
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
                        cutTorque((min(MUFiring2{i})*2048-auxTimeCut(1))),cutTorque((max(MUFiring2{i})*2048-auxTimeCut(1))),...
                        cutTorque((min(MUFiring2{j})*2048-auxTimeCut(1))),cutTorque((max(MUFiring2{j})*2048-auxTimeCut(1))),...
                        cutEMG((min(MUFiring2{i})*2048-auxTimeCut(1))),cutEMG((max(MUFiring2{i})*2048-auxTimeCut(1))),...
                        cutEMG(min(MUFiring2{j})*2048-auxTimeCut(1)),cutEMG((max(MUFiring2{j})*2048-auxTimeCut(1)))];%%%%%%%%%%
                    c=c+1;
                end
                
                if keepPics >= 3
                    %savename = [cleanedFile(1:end-4),'_Unit_', num2str(TRUEMU(i), '%02i'),'_Unit_',num2str(TRUEMU(j), '%02i'),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))];
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
        for i = 1:length(MUFiring2)
            starttime = MUFiring2{i}(1);
            endtime = MUFiring2{i}(end);
            f99 = figure(99);clf;
            
            if c<21
                plot(MUFiring2{i}(2:end),1./(diff(MUFiring2{i})),'o','Color',Cmap(c,:))
                hold on;
                plot(starttime:1/fsamp:endtime,MUTrain(i,starttime*fsamp:endtime*fsamp),'Color',Cmap(c,:),'LineWidth', 2)
            else
                c=1;
                plot(MUFiring2{i}(2:end),1./(diff(MUFiring2{i})),'o','Color',Cmap(c,:))
                hold on;
                plot(starttime:1/fsamp:endtime,MUTrain(i,starttime*fsamp:endtime*fsamp),'Color',Cmap(c,:),'LineWidth', 2)
            end
            title({SubjectID, ['Unit ', num2str(TRUEMU(i), '%02i')]}, 'Interpreter', 'none');
            axis(plot_axis)
            
            if keepPics >= 2
                savenamef99 = [SubjectID,'_Unit_', num2str(TRUEMU(i), '%02i'),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))];
                print(f99,'-dpdf',fullfile(pdfdir, savenamef99))
            end
            c=c+1;
        end
    end
    
    %%% uncertain of purpose of this section as of 1/29/22
    [M,I] = max(cutTorque(5*fsamp:end));
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
    Torqueout(10) = max(cutEMG);
    
    output(size(output,1)+1,:) = Torqueout;
    
    %%% plot and save sumary figure
    if Pics >= 1
        try
            close(sumFig)
        catch
        end
        sumFig = figure(101); hold all
        pos = get(gcf,'position');
        set(gcf,'position',[pos(1) pos(2)-pos(4) pos(3) pos(4)*2])
        subplot(8,2,[1 4]); hold all
        xlim([MUTime(1) MUTime(2)]);
        ylim([-10 20]);
        
        c = 1;
        
        for i = 1:length(windowDR)
            try
                plot(windowST{i},windowDR{i},'LineWidth',1.5,'Color',Cmap(c,:));
                c=c+1;
            catch
                c=1;
                plot(windowST{i},windowDR{i},'LineWidth',1.5,'Color',Cmap(c,:));
            end
            title({[clusterFile(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],['All ',num2str(length(MUFiring2)),' units']},'interpreter','none')
        end
        
        try
            plot(timeTorque,dfTQfeedback,'k','LineWidth',2)%
            plot(timeTorque,pfTQfeedback,'r','LineWidth',2)%
            plot(timeTorque,newEMGFeedback,'SeriesIndex', agonistEMGColor,'LineWidth',1,'LineStyle',':')
        catch
            plot(timeTorque,cutTorque,'Color','k','LineWidth',2)
            plot(timeTorque,cutEMG-10,'Color', [.3 .3 .3],'LineWidth',1.5,'LineStyle',':')
        end
        plot([1/fsamp:1/fsamp:size(MUTrain,2)/fsamp],MUTrain(size(MUTrain,1),:),'Color', [.1 .1 .1],'LineWidth',1.5)
        
        title4All = [SubjectID,'_allunits_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp',num2str(jj),'.pdf'];
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
            sumFigTitle = [SubjectID, '_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp', num2str(jj),'.pdf'];
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
        
        dlmwrite(fullfile(exceldir,[SubjectID,'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_Ramp', num2str(jj),'_deltaf.xls']),output,'delimiter','\t')
        
    end
end
end

