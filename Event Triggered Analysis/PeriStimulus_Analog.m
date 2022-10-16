function [periAnalogTime,periAnalogSignal] = PeriStimulus_Analog(AnalogSignal, triggerTimeStamps, timeAnalog)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
loopAnalog = [];
    for i = 1:length(triggerTimeStamps)
        windowAnalog = AnalogSignal(timeAnalog>triggerTimeStamps(i)-2 & timeAnalog<triggerTimeStamps(i)+2);
        loopAnalog = [loopAnalog; windowAnalog];
    end
    
    periAnalogSignal = loopAnalog;
    periAnalogTime = linspace(-1,1.5,size(periAnalogSignal,2));
end

