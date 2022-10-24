
function [allbadchan, data] = Preprocess_OTB_EMG(arrayMatrix, Muscle, ArrayChannels, TotalChans, fsamp) 

%% DESCRIPTION

%   Loads OT Biolab .sig file.  Plots total number of inputs;
%   multi-channel matrices are plotted in groups of 16. Prompts user for
%   bad channels to remove.  Outputs the bad channels.

% UPDATED   Matt Topley
% AUTHOR    Laura Miller
% BUTCHER   Christopher Thompson

% DATE CREATED      6-Feb-2013
% DATE BUTCHERED   23-Mar-2013
% DATE Updated   08-Oct-2022

%% Read sig file - This is different depending on version of OTB Lab
x_coordstart = 1;
allbadchan = [];

%     %%%% opening the sig file is different with OTB+ files. The script is updated to open old and new OTB files %%%%
% [~, Trialname] = fileparts(TrialPath);
% 
% sigFile = fullfile(TrialPath, [Trialname,'.sig']);
% 
% if ~exist(sigFile)
%     % this will only work with the new OTB+ files
% sigFile = dir(fullfile(TrialPath,'*.sig'));
% sigFile = fullfile(sigFile.folder, sigFile.name);
% end 
% 
% f = fopen(sigFile);
% 
% %%% extract data from Quattro 
% data = fread(f,[TotalChans+8,inf],'short');
data = double(arrayMatrix(:,x_coordstart:end));

%% plotting EMG for removal of bad channels
emgFig = figure(101);
for k = 1:length(ArrayChannels)
    current  = ArrayChannels{k};
    if length(current) > 1
        
        %Array channels broken into groups of 16 
        current16 = {[current(1:16)] [current(17:32)] [current(33:48)] [current(49:64)]};

        for m = 1:length(current16)
            
                channelNumber = string(current16{m});
                
                for i = 1:16
                    plotdata(current16{m}(i),:) = data(current16{m}(i),:) - ones(1,length(data))*(i-1)*1000;
                end
                plot((1:length(plotdata))/fsamp,plotdata(current16{m},:)');
                
                title([Muscle, ' - ', num2str(current16{m}(i-15)), ' : ', num2str(current16{m}(i))]);
                yticks([-15000:1000:1000]);
                yticklabels([flip(channelNumber), ""]);
                emgYAxis = get(gca,'ylim');
                ylim([emgYAxis(1) 1000])
                ylabel('Channel Number');
                
                %% Remove bad channels
                badchan = input('Enter bad channels or empty matrix: e.g. [1,2]: ');
                chansEntered = ismember(badchan, current16{m});
                
                % check if the bad channels need to be reentered
                if isempty(badchan)
                    
                elseif max(chansEntered) == 0 
                badchan = input('Re-enter bad channels or empty matrix: e.g. [1,2]: ');
                end 
                
                if badchan == 'n'
                    allbadchan = badchan;
                else
                    allbadchan = [allbadchan; badchan'];
                end
                clf(emgFig); % Clear Figure
        end
        
    else
        plotdata(current,:) = data(current,:);
        if badchan == 'n'
            allbadchan = badchan;
        else
            allbadchan = [allbadchan; badchan'];
        end
        clf
    end
end

%%%%%%%%% Differential recordings, comment out if Monopolar
if max(ArrayChannels{1}) == 64 && min(ArrayChannels{1}) == 1
    allbadchan = [allbadchan; 64];
elseif max(ArrayChannels{1}) == 128 && min(ArrayChannels{1}) == 65
    allbadchan = [allbadchan; 128];
elseif max(ArrayChannels{1}) == 192 && min(ArrayChannels{1}) == 129
    allbadchan = [allbadchan; 192];
elseif max(ArrayChannels{1}) == 256 && min(ArrayChannels{1}) == 193
    allbadchan = [allbadchan; 256];
elseif max(ArrayChannels{1}) == 320 && min(ArrayChannels{1}) == 257
    allbadchan = [allbadchan; 320];
elseif max(ArrayChannels{1}) == 384 && min(ArrayChannels{1}) == 321
    allbadchan = [allbadchan; 384];
end

allbadchan = sort(allbadchan);

%fclose all;
end
