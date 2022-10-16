%%
%  clear all;
clear vars 

trial = 3;
contraction = 'EMG';
filenames = dir('*CoCo_coher*.mat');

% Tableaumap = Tableau;
for i = 1:length(filenames)
    try
        load(filenames(i).name,'taF', 'taz', 'taCOFz', 'taCOF2z', 'taCOF','MULength', 'pF', 'pz', 'pCOFz', 'pCOF2z')
        
        
        deltaband = trapz(taF(1:50), taz(1:50));    % 0 - 5Hz
        alphaband = trapz(taF(50:150), taz(50:150));    % 5 - 15Hz
        betaband = trapz(taF(150:350), taz(150:350));    % 15 - 35Hz
        
        
        
        legend('AutoUpdate', 'off')
        hold on
        area(taF(1:50), taz(1:50), 'SeriesIndex', 2, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        area(taF(50:150), taz(50:150), 'SeriesIndex', 4,'FaceAlpha', 0.3, 'EdgeColor', 'none')
        area(taF(150:350), taz(150:350),'SeriesIndex', 6, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        c1 = plot(taF,ones(1,length(taF))*taCOFz,'r');
        c2 = plot(taF,ones(1,length(taF))*taCOF2z,'k');
        hold off
        legend('-DynamicLegend', 'Location', 'northeast', 'AutoUpdate', 'on', 'FontSize', 15)
        bands = sprintf(' Delta: %0.2f Alpha: %0.2f Beta: %0.2f',deltaband,alphaband,betaband);
        %trialname = [contraction, ' Delta: ', num2str(deltaband), ', Alpha: ', num2str(alphaband), ', Beta: ', num2str(betaband)];
        trialname = [contraction, bands]
        hold on
        plot(taF,taz,'SeriesIndex', trial, 'LineWidth', 2, 'DisplayName', trialname)
        hold off
        xlim([0 50])
        ylim([0 10])
        t =gca;
        t.XAxis.FontSize  = 15;
        t.YAxis.FontSize  = 15;
        title('TA-TA');
        
        t.XLabel.String = 'Frequency (Hz)';
        t.YLabel.String = 'Coherence (z-transform)';
    end


% if trial == 1
% title('No Coco')
% elseif trial == 2
%     title('With Coco')
% % else
% %     title('CoCo hold TA-MG Coherence')
% end 

  print (gcf,'-dpdf', '-bestfit',[filenames(i).name(1:end-20), contraction, '_TA_Coherence.pdf']);
end
%%  MG single subject plots
%%
%  clear all;
clear vars 

trial = 3;
contraction = 'EMG';
filenames = dir('*CoCo_coher*.mat');

% Tableaumap = Tableau;
for i = 1:length(filenames)
   try
%     load(filenames(i).name,'taF', 'taz', 'taCOFz', 'taCOF2z', 'taCOF','MULength', 'pF', 'pz', 'pCOFz', 'pCOF2z')
        load(filenames(i).name,'mgF', 'mgz', 'mgaCOFz', 'mgCOF2z', 'mgCOF','MULength', 'pF', 'pz', 'pCOFz', 'pCOF2z')   
   end 
end
 count = i;
 try
 deltaband = trapz(mgF(1:50), mgz(1:50));    % 0 - 5Hz
 alphaband = trapz(mgF(50:150), mgz(50:150));    % 5 - 15Hz
 betaband = trapz(mgF(150:350), mgz(150:350));    % 15 - 35Hz
 
 
    
legend('AutoUpdate', 'off')
hold on
area(mgF(1:50), mgz(1:50), 'SeriesIndex', 2, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
area(mgF(50:150), mgz(50:150), 'SeriesIndex', 4,'FaceAlpha', 0.3, 'EdgeColor', 'none')
area(mgF(150:350), mgz(150:350),'SeriesIndex', 6, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
c1 = plot(mgF,ones(1,length(mgF))*mgCOFz,'r');
c2 = plot(mgF,ones(1,length(mgF))*mgCOF2z,'k');
hold off
legend('-DynamicLegend', 'Location', 'northeast', 'AutoUpdate', 'on', 'FontSize', 15)
bands = sprintf(' Delta: %0.2f Alpha: %0.2f Beta: %0.2f',deltaband,alphaband,betaband);
 %trialname = [contraction, ' Delta: ', num2str(deltaband), ', Alpha: ', num2str(alphaband), ', Beta: ', num2str(betaband)]; 
 trialname = [contraction, bands]
 hold on
    plot(mgF,mgz,'SeriesIndex', trial, 'LineWidth', 2, 'DisplayName', trialname)
    hold off
xlim([0 50])
ylim([0 10])
t =gca;
t.XAxis.FontSize  = 15;
t.YAxis.FontSize  = 15;
title('MG-MG');

t.XLabel.String = 'Frequency (Hz)';
t.YLabel.String = 'Coherence (z-transform)';
 end 

% if trial == 1
% title('No Coco')
% elseif trial == 2
%     title('With Coco')
% % else
% %     title('CoCo hold TA-MG Coherence')
% end 

 print (gcf,'-dpdf', '-bestfit',[contraction, '_MG_Coherence.pdf']);