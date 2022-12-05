function [goodUnits, TRUEMU, rSpikes] = Remove_MUs_User(MotorUnits)
%%% This function takes the file, the muscles motor units, the smoothed binary spike
%%% train, and the contraction window to help remove bad motor units

% the function will calculate the COV of the smoothed spike trains within the contraction window
% and plot them in the legend. Two windows will plot the units and spike trains for better selection
% bad units will be removed and one of the windows will be closed. The
% other window will stay open and be reused later in main script
varC =  class(MUPulses);
fsamp = 2048;
figure

if varC == 'cell'
    for i = 1:length(MotorUnits)
        hold on
        plot(MotorUnits{i}(2:end) ./fsamp, 1./diff(MotorUnits{i}) .* fsamp, 'o')
        hold off
    end
    
else
    if size(MotorUnits, 2) > size(MotorUnits, 1)
        plot(MotorUnits')
    else
        plot(MotorUnits)
    end
end

rSpikes = input('Enter the spikes to remove');

goodUnits = MotorUnits;

TRUEMU = 1:length(MotorUnits);
goodUnits(rSpikes) = [];
TRUEMU(rSpikes) = [];
fig = gcf;
close(fig)
end
