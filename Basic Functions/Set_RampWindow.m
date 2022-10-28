function [rampStart] = Set_RampWindow(varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048; % always declare this at the beginning of our functions to be safe
 
condition = varargin{1};
if varargin >= 2
    trace = varargin{2}; % parsing inputs
end

if varargin >= 3
    Torque = varargin{3};  
end

if varargin == 4
    EMGFeedback = varargin{4};
end

if ~contains(condition, 'coco') && ~ischar(trace)
    trace = TraceFeedback-mean(TraceFeedback(1:fsamp));
    trace = movmean(trace,fsamp);
    risingEdge = trace>2;
    [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 50000);
    rampStart = round(rampStart./fsamp)-5;
    
elseif ~contains(condition, 'coco') && ischar(trace)
     cop = Torque-mean(Torque(1:fsamp));
    cop = movmean(cop,fsamp);
    risingEdge = cop>2;
    [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 50000);
    rampStart = round(rampStart./fsamp)-5;
    
elseif ~contains(condition, 'coco') && ischar(trace)  && ischar(Torque)
     findEMGRamp = movmean(EMGFeedback, fsamp);
    risingEdge = findEMGRamp>2;
    [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 50000);
    if length(rampStart) < 3
        [~,rampStart] = findpeaks(risingEdge*2, 'MinPeakProminence', 1,'MinPeakDistance', 25000);
    end
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

