function [MotorUnits] = SortUnits(MotorUnits)
%This function sorts motor units from high threshold to small threshold, or
%shortest spike trains to longest

%   "Empty" units are automatically removed as well
MotorUnits(find(cellfun(@isempty,MotorUnits))) = [];

[~,FiringLen] = sort(cellfun(@length,MotorUnits));
MotorUnits = MotorUnits(FiringLen);

end

