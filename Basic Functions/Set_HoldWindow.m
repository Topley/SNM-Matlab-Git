function [holdStart] = Set_HoldWindow(TraceFeedback,AnalogFeedback)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;

if mean(TraceFeedback) > 2
risingEdge = TraceFeedback>max(TraceFeedback)*.75;
[~,holdStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'NPeaks', 1);
holdStart = round(holdStart./fsamp)+5;

else
%findEMGRamp = movmean(EMGFeedback, fsamp);
risingEdge = AnalogFeedback>max(AnalogFeedback)*.75;
[~,holdStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'NPeaks', 1);
holdStart = round(holdStart./fsamp)+5;


end

