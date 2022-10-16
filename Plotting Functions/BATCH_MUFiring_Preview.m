
function BATCH_plotMUFiring_AllUnits
filenames = dir('*TA_5_v23decomposed.mat')
% filenames = dir('*Sol_5_v23decomposed_MUCLEANED.mat')

%filenames = ls('CoCoTest09_44426_RTA_1_v23decomposed_MUCLEANED.mat')
for i = 1:size(filenames,1)
    try
        plot_AllMUFiring(filenames(i).name)
    end
    close all
end
end


function [] = plot_AllMUFiring(filenameTA)
% plotMUFiring    loads motor units decomposed from cluster; plots units
%
% USAGE
% check = 0;
% if isfile(strcat(filename(1:end-4), '_CHECK.pdf')) == 1
%     return
% end
% checkFile = replace(filenameTA, 'Sol', 'TA');
% if exist(checkFile) == 2
%     return
%     %disp('Yesy')
% end

fsamp = 2048;

[tableau] = Tableau;

  clusterFileVariables = {'fsamp','xCorChannel'};
% clusterFileVariables = {'fsamp','CrossCorrEMG'};
% filenameTA = replace(filenameTA, 'Sol', 'TA');
try
    TAPulses = load(filenameTA, 'MUPulses');
    TAPulses = TAPulses.MUPulses;
    [~,DRlen] = sort(cellfun(@length,TAPulses));
    TAFiring = TAPulses(DRlen);
end

taUnder = strfind(filenameTA,'_');
if ~contains(filenameTA,'MUCLEANED')
    taClusterFile = [filenameTA(1:(taUnder(end)-1)),'.mat'];
else
    taClusterFile = [filenameTA(1:(taUnder(end-1)-1)),'.mat'];
end

try
    taStruct =  load(taClusterFile, clusterFileVariables{:});
end

filenameSol = replace(filenameTA,'TA','Sol');
filenameMG = replace(filenameTA,'TA','MG');
% filenameLG = replace(filenameTA,'TA','LG');

try
    SolPulses = load(filenameSol, 'MUPulses');
    SolPulses = SolPulses.MUPulses;
    [~,DRlen] = sort(cellfun(@length,SolPulses));
    SolFiring = SolPulses(DRlen);
end

try
    MGPulses = load(filenameMG, 'MUPulses');
    MGPulses = MGPulses.MUPulses;
    [~,DRlen] = sort(cellfun(@length,MGPulses));
    MGFiring = MGPulses(DRlen);
end

% try
%     LGPulses = load(filenameLG, 'MUPulses');
%     LGPulses = LGPulses.MUPulses;
%     [~,DRlen] = sort(cellfun(@length,LGPulses));
%     LGFiring = LGPulses(DRlen);
% end

try
    %%NORMAL
    Torque = taStruct.xCorChannel;
catch
    Torque = taStruct.CrossCorrEMG;
end
Torque = (Torque-mean(Torque(1:fsamp)))./100;
[B,A] = butter(3,[10]*2/fsamp,'low');
Torque_butter = filtfilt(B,A,Torque')';
Torque_time = 1/fsamp:1/fsamp:length(Torque)/fsamp;

h=figure;
unitIdx = 0;
try
    for j = 1:size(TAFiring,2)
        TAFiring_loop = [];
        TAFiring_loop = TAFiring{j};
        unitIdx = unitIdx+1;
        TAISI = diff(TAFiring_loop)/fsamp; % Interspike interval in sec
        TAIDR = 1./TAISI; % instantaneous discharge rate
        TATime = TAFiring_loop(2:end)/fsamp; % Time of discharge in sec
        
        subplot(2,1,2);
        plot(TATime,TAIDR-(unitIdx-1)*20,'.','MarkerSize',8,'SeriesIndex',1);hold all    % offset axis
        xlim([0 max(Torque_time)]);
    end
end

try
    for ij = 1:size(SolFiring,2)
        unitIdx = unitIdx+1;
        SolFiring_loop = [];
        SolFiring_loop = SolFiring{ij};
        
        SolISI = diff(SolFiring_loop)/fsamp; % Interspike interval in sec
        SolIDR = 1./SolISI; % instantaneous discharge rate
        SolTime = SolFiring_loop(2:end)/fsamp; % Time of discharge in sec
        
        subplot(2,1,2);
        plot(SolTime,SolIDR-(unitIdx-1)*20,'.','MarkerSize',8,'SeriesIndex',7);hold all    % offset axis
        xlim([0 max(Torque_time)]);
    end
end

try
    for jj = 1:size(MGFiring,2)
        MGFiring_loop = [];
        MGFiring_loop = MGFiring{jj};
        unitIdx = unitIdx+1;
        MGISI = diff(MGFiring_loop)/fsamp; % Interspike interval in sec
        MGIDR = 1./MGISI; % instantaneous discharge rate
        MGTime = MGFiring_loop(2:end)/fsamp; % Time of discharge in sec
        
        subplot(2,1,2);
        plot(MGTime,MGIDR-(unitIdx-1)*20,'.','MarkerSize',8,'SeriesIndex',5);hold all    % offset axis
        xlim([0 max(Torque_time)]);
    end
end

try
    subplot(2,1,1)
    plot(Torque_time,Torque);hold all
    plot(Torque_time,Torque_butter,'LineWidth',2);
    axis([0 max(Torque_time) min(Torque) max(Torque)]);
end

p = get(subplot(2,1,1),'position');
p(2) = p(2)*1.4;  % bottom
p(4) = p(4)*0.4;  % height
set(subplot(2,1,1), 'position', p);
try
    xlim([0 max(Torque_time)]);
end
box off

title([filenameTA(1:(taUnder(2))),'AllUnits_v23decomposed.mat'])
%title(filenameTA(1:end-4),'interpreter','none')

subplot(2,1,2);

y = get(subplot(2,1,2),'ylim');

if y(2)>100
    y(2)=100;
end


% if ~isempty(TAFiring)
    
    y(1) = -(unitIdx-1).*20;
    
    set(subplot(2,1,2),'ylim',y)
    
    p = get(subplot(2,1,2),'position');
    p(2) = p(2)*.1; % Subtracts 90 percent from bottom
    p(4) = p(4)*2.3;  % Add 200 percent to height
    set(subplot(2,1,2), 'position', p);
    set(gca,'YGrid','on');
    set(gca,'GridLineStyle','-');
    set(gca, 'YColor', [.4, .4, .4]);
    
    Ticklines = [unitIdx*20*-1:20:20];
    tickLabels = [1:length(Ticklines)]';
    
    set(gca,'YTick',Ticklines, 'YTickLabels',tickLabels);
% end
box off

savefilename = strcat([filenameTA(1:(taUnder(2))),'AllUnits_v23decomposed'], '_CHECK.pdf');
print (h,'-dpdf',savefilename);
pause(1)
close(h)
%
% if check == 1
%     if isfile(strcat(filenameTA(1:end-4), '_MUCLEANED_REVIEW.pdf')) == 1;
%         close(h)
%         return
%     end
%     savefilename = strcat(filenameTA(1:end-4), '_CHECK.pdf');
%     print (h,'-dpdf',savefilename);
%     pause(1)
%     close(h)
%
% elseif check == 0
%     try
%         warning('off', 'MATLAB:DELETE:FileNotFound');
%         delete(strcat(filenameTA(1:end-4),'CHECK.pdf'));
%     catch
%     end
%     %      if isfile(strcat(filename(1:end-4), '_REVIEW.pdf')) == 1;
%     %         close(h)
%     %         return
%     %     end
%     savefilename = strcat(filenameTA(1:end-4), '_REVIEW.pdf');
%     print (h,'-dpdf',savefilename);
%     pause(1)
%     close(h)
% else
% end

end

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
    0.6196    0.8549    0.8980]; % load colormap
end


