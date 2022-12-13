function [segmentRotMat] = RotationMatrix(segment, markers)
%%%% Rotation Matrix 
%   First argument specifies the segment the rotation matrix represents.
%   The second argument is a cell array of the markers to build the matrix
%   from

if contains(segment, 'Hip') && length(markers) < 3 
markers{3} = (markers{2} - markers{1}) ./ 2;
end 

marker1 = markers{1};
marker2 = markers{1};
marker3 = markers{3};

try
    apOffsetVect = marker3 - marker1;
    
    axisZ = marker2 - marker1;
    axisY = cross(apOffsetVect, axisZ);
    axisX = cross(axisZ, axisY);
    
    Xvect = normr(axisX);
    Yvect = normr(axisY);
    Zvect = normr(axisZ);
    
    frames = size(Xvect,1);
    xRot = reshape(Xvect,1, 3, frames);
    RotMat = cat(1,xRot, [reshape(Yvect,1, 3, frames);reshape(Zvect,1, 3, frames)]);
    segmentRotMat = permute(RotMat,[2 1 3]);
catch
    disp(['Error creating', segment, 'Segment Rotation Matrix'])
end

end

