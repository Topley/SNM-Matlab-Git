function [unitsISI] = ISIs(MUFiring)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
unitISI = cell(1,length(MUFiring));

for i = 1:length(MUFiring)
    unitsISI{i} = diff(MUFiring{i})/fsamp;
end

end

