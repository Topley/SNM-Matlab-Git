function [smoothIDRs,smoothCST,smoothTime] = Pulses2Firing(MotorUnits, AnalogReference)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%% This is where the units are filtered
fsamp = 2048;
for i = 1:length(MotorUnits)
    MuF = MotorUnits{i};
    %index = round(MuF*fsamp);
    firing(i,MuF) = 1;
end

firing(:,length(firing)+fsamp)=0; %adds a second of 0s to the end
smoothIDRs = 1./sum(hanning(fsamp*2))*filtfilt(hanning(fsamp*2),1,firing')';%AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
smoothCST = mean(smoothIDRs);
smoothTime = [0:length(firing)-1]/fsamp;

end

