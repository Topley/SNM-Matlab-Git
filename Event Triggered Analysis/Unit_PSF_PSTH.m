function [periStimUnitTime,periStimUnitIDR] = Unit_PSF_PSTH(MotorUnit,Flicks, stimulusWindow)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
windowStart = stimulusWindow(1);
windowEnd = stimulusWindow(2);

totalUnitLoops = [];
normMUFiring = 1./diff(MotorUnit./fsamp);
normUnitTime = MotorUnit(2:end)./fsamp;
for i = 1:length(Flicks)
    loopUnit = [];
    loopUnit(2,:) = normUnitTime(normUnitTime>Flicks(i)+windowStart & normUnitTime<Flicks(i)+windowEnd) - Flicks(i);
    loopUnit(1,:) = normMUFiring((normUnitTime>Flicks(i)+windowStart).*fsamp & (normUnitTime<Flicks(i)+windowEnd).*fsamp);
    totalUnitLoops = [loopUnit totalUnitLoops];
end

periStimUnitTime = totalUnitLoops(2,:);
periStimUnitIDR = totalUnitLoops(1,:);

end

