function [normFiring] = normMUFiring(MUFiring)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
if size(MUFiring,2) > size(MUFiring,1)
    MUFiring = MUFiring';
end 

[b,a] = butter(4, 1/round(fsamp/2.5),'low');
filtMUFiring = filtfilt(b,a,MUFiring);

for i = 1:size(filtMUFiring,2)
    for j = 1:length(filtMUFiring(:,i))
        normFiring(j,i) = (filtMUFiring(j) - min(filtMUFiring(:,i))) / (max(filtMUFiring(:,i)) - min(filtMUFiring(:,i)));
    end
end

end

