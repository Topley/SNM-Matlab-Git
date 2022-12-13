function [segmentRotMat] = TrunkRotationMatrix(Markers)
%%%% Rotation Matrix
%   First argument specifies the segment the rotation matrix represents.
%   The second argument is a cell array of the markers to build the matrix
%   Markers shouold be entered in order of Lateral first, Medial second,
%   and the offset marker third.

marker1 = Markers.STR;
marker2 = Markers.C7;

try
    marker3 = Markers.SAC;
catch
    marker3 = (Markers.RPSI - Markers.LPSI) ./ 2;
end

axisY = marker2 - marker3;
apOffsetVect = marker1 - marker2;
axisZ = cross(apOffsetVect, axisY);
axisX = cross(axisY, axisZ);

Xvect = normr(axisX);
Yvect = normr(axisY);
Zvect = normr(axisZ);

frames = size(Xvect,1);
xRot = reshape(Xvect,1, 3, frames);
RotMat = cat(1,xRot, [reshape(Yvect,1, 3, frames);reshape(Zvect,1, 3, frames)]);
segmentRotMat = permute(RotMat,[2 1 3]);

end

