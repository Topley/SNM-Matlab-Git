function [axisRange, checkFig] = Set_RampWindow_Click(filename,checkFig, SpikeTrain, window)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 % Click analysis range for Antagonist
 fsamp =2048;
 
 try
        checkFig = gcf;
        clf('reset')
        hold on
        fill([window(1), window(2),window(2),window(1)],[0 0 25 25],'k','facealpha',0.1,'edgecolor','none');
        plot([0:length(SpikeTrain)-1]/fsamp,2*fftfilt(hanning(fsamp),SpikeTrain'));
        hold off
        title(filename);
        disp('Choose analysis window ')
        [ax, ~] = ginput(2);
        axisRange = round(ax);
        
 catch
        disp(['Error in analysisWindowClick, plotting ', filename, ' MUs']);
 end
%  close(checkFig)
end

