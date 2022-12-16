function [RotMat, ASISBreadth, HipOrigin, JointCenter] = HipRotationMatrix(Markers)
%%%% Rotation Matrix
%   First argument specifies the segment the rotation matrix represents.
%   The second argument is a cell array of the markers to build the matrix
%   Markers shouold be entered in order of Lateral first, Medial second,
%   and the offset marker third.

TribalSys = [0,1,0; 1,0,0;0,0,1];
rHipConstants = [-24; 32; -21]; % as a percentage
lHipConstants = [-24; -32; -21]; % as a percentage

marker1 = Markers.LASI;
marker2 = Markers.RASI;

try
    marker3 = Markers.SAC;
catch
    marker3 = (Markers.RPSI - Markers.LPSI) ./ 2;
end

HipOrigin = (marker2 + marker1) ./2;

apOffsetVect = marker3 - marker1;
axisZ = marker2 - marker1;
axisY = cross(apOffsetVect, axisZ, 2);
axisX = cross(axisY,axisZ, 2);

Xvect = normr(axisX);
Yvect = normr(axisY);
Zvect = normr(axisZ);

ASISBreadth = mean(Zvect);

xRot = reshape(Xvect, 3,1, []);
yRot = reshape(Yvect, 3,1, []);
zRot = reshape(Zvect, 3,1, []);

RotMatrix = [xRot,yRot,zRot];
HipRotMat = pagemtimes(RotMatrix, TribalSys);

%HipRotMat = permute(RotMatrix,[2 1 3]);
%HipRotMat = RotMatrix;
RotMat.Pelvis = HipRotMat;

locateRight = ASISBreadth .* (rHipConstants ./ 100)';
rotateRight = pagemtimes(HipRotMat, 'transpose', locateRight, 'transpose');
transposeRight = permute(rotateRight, [3 1 2]);

locateLeft = ASISBreadth .* (lHipConstants ./ 100)';
rotateLeft = pagemtimes(HipRotMat, 'transpose', locateLeft, 'transpose');
transposeLeft = permute(rotateLeft, [3 1 2]);

JointCenter.RHip = transposeRight + HipOrigin; 
JointCenter.LHip = transposeLeft + HipOrigin; 

end

