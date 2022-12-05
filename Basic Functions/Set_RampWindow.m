function [rampStart, rampEnd] = Set_RampWindow(varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048; % always declare this at the beginning of our functions to be safe

AuxSignal = varargin{1};
baseAuxSignal = AuxSignal-mean(AuxSignal(1:fsamp));
auxSignal = movmean(baseAuxSignal, 2*2048);
thresh = 2;
ramps = diff(auxSignal > thresh);
trialTime = [0:length(auxSignal)-1]./fsamp;
rampStart1 = round(trialTime(ramps == 1) - 2);
rampEnd1 = trialTime(ramps == -1) + 2;

for i = 1:length(rampStart1)
    [~, idx] = max(auxSignal(fsamp*rampStart1(i) : fsamp*rampEnd1(i)));
    upR = round(idx/fsamp - 10);
    if upR > 1
        rampStart(i) = upR;
        rampEnd(i) = round(idx/fsamp + 10);
    else
        rampStart(i) = 1;
        rampEnd(i) = 21;
    end
end
if nargin > 1
    AuxSignal2 = varargin{2};
    baseAuxSignal2 = AuxSignal2-mean(AuxSignal2(1:fsamp));
    auxSignal2 = movmean(baseAuxSignal2, 2*2048);
    thresh2 = 2;
    ramps2 = diff(auxSignal2 > thresh2);
    trialTime2 = [0:length(auxSignal2)-1]./fsamp;
    rampStart2 = round(trialTime2(ramps2 == 1) - 2);
    rampEnd2 = rampStart2 + 30;
    %rampEnd2 = trialTime2(ramps2 == -1);
end
%
% rampWidth1 = abs(mean(rampStart1));
% rampWidth2 = abs(mean(rampStart2));
%     if rampWidth2 < rampWidth1
%         rampStart = rampStart2;
%         rampEnd = rampEnd2;
%     else
% rampStart = rampStart1;
% rampEnd = rampEnd1;
%     end

end

