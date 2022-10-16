function [unitsISI] = ISIs(MUFiring)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    fsamp = 2048;
    unitsISI = diff(MUFiring)/fsamp;
end

