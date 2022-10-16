function [CST_Time, CST_IDR, perStimulusCST_Time, perStimulusCST_IDR] = PeriStimulus_CST(MotorUnits, Flicks, stimulusWindow)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
windowStart = stimulusWindow(1);
windowEnd = stimulusWindow(2);

CST = sort(cell2mat(MotorUnits));
CST_Time = [];
CST_IDR = [];
for i = 1:length(MotorUnits)
    loopMU = [];
    loopMU = (1./diff((MotorUnits{i}./fsamp)));
    
    CST_IDR = [loopMU, CST_IDR];
    CST_Time =  [(MotorUnits{i}(2:end)./fsamp), CST_Time];
end

loopCSTPSF = [];
for i = 1:length(Flicks)
    loopPSF = [];
    loopPSF(1,:) = CST_IDR(CST_Time>Flicks(i)+windowStart & CST_Time<Flicks(i)+windowEnd);
    loopPSF(2,:) = CST_Time(CST_Time>Flicks(i)+windowStart & CST_Time<Flicks(i)+windowEnd) - Flicks(i);
    loopCSTPSF = [loopPSF loopCSTPSF];
end

perStimulusCST_Time = loopCSTPSF(2,:);
perStimulusCST_IDR = loopCSTPSF(1,:);



end

