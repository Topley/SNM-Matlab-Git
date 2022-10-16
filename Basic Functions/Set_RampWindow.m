function [rampStart] = Set_RampWindow(TraceFeedback,EMGFeedback, EMG, condition)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
if condition == 0
    trace = TraceFeedback-mean(TraceFeedback(1:fsamp));
    trace = movmean(trace,fsamp);
    risingEdge = trace>2;
    [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 50000);
    rampStart = round(rampStart./fsamp)-5;
else
    findEMGRamp = movmean(EMGFeedback, fsamp);
    risingEdge = findEMGRamp>2;
    [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 50000);
    if length(rampStart) < 3
        [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 25000);
    end
    rampStart = round(rampStart./fsamp)-5;
    
end

end

