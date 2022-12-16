function [segRotMat, JointCenter] = segmentRotationMatrix(segment, Markers, JointCenter)
%%%% Rotation Matrix
%   First argument specifies the segment the rotation matrix represents.
%   The second argument is a cell array of the markers to build the matrix
%   from
TribalSys = [0,1,0; 1,0,0;0,0,1];
segmentName = segment;
    
switch segmentName
    case 'RThigh'
        latMarker = Markers.RKNE;
        medMarker = Markers.RMK;
        ProximalJC = JointCenter.RHip;
        try
        DistalJC = JointCenter.RKNE;
        catch
        DistalJC = (latMarker + medMarker) ./2;
        JointCenter.RKNE = DistalJC;
        end 
        
    case 'LThigh'
        latMarker = Markers.LKNE;
        medMarker = Markers.LMK;
        ProximalJC = JointCenter.LHip;
        try
        DistalJC = JointCenter.LKNE;
        catch
        DistalJC = (latMarker + medMarker) ./2;
        JointCenter.LKNE = DistalJC;
        end 
        
    case 'RShank'
        latMarker = Markers.RANK;
        medMarker = Markers.RMA;
        ProximalJC = JointCenter.RKNE;
        try
        DistalJC = JointCenter.RANK;
        catch
        DistalJC = (latMarker + medMarker) ./2;
        JointCenter.RANK = DistalJC;
        end 
        
    case 'LShank'
        latMarker = Markers.LANK;
        medMarker = Markers.LMA;
        ProximalJC = JointCenter.LKNE;
        try
        DistalJC = JointCenter.LANK;
        catch
        DistalJC = (latMarker + medMarker) ./2;
        JointCenter.LANK = DistalJC;
        end 
        
    case 'RFoot'
        toeMarker = Markers.RTOE;
        heelMarker = Markers.RHEE;
        ankleJC = JointCenter.RANK;
        
    case 'LFoot'
        toeMarker = Markers.LTOE;
        heelMarker = Markers.LHEE;
        ankleJC = JointCenter.LANK;
        
    otherwise
        disp('No segment chosen')
end

if ~contains(segmentName, 'Foot')
    axisY = ProximalJC - DistalJC;
    tempVect = latMarker - medMarker;
    axisX = cross(axisY, tempVect, 2);
    axisZ = cross(axisX, axisY, 2);
    
    Xvect = normr(axisX);
    Yvect = normr(axisY);
    Zvect = normr(axisZ);
    
    xRot = reshape(Xvect,3, 1, []);
    yRot = reshape(Yvect,3, 1, []);
    zRot = reshape(Zvect,3, 1, []);
    RotMat = [xRot,yRot,zRot];

    segRotMat = pagemtimes(RotMat, TribalSys);
    %segRotMat = RotMat;
    %RotMat = cat(1,xRot, [reshape(Yvect,1, 3, frames);reshape(Zvect,1, 3, frames)]);
    %segRotMat = permute(RotMat,[2 1 3]);
else
    axisX = toeMarker - heelMarker;
    tempVect = ankleJC - heelMarker;
    axisZ = cross(axisX, tempVect);
    axisY = cross(axisZ, axisX);
    
    Xvect = normr(axisX);
    Yvect = normr(axisY);
    Zvect = normr(axisZ);
    
    xRot = reshape(Xvect,3, 1, []);
    yRot = reshape(Yvect,3, 1, []);
    zRot = reshape(Zvect,3, 1, []);
    RotMat = [xRot,yRot,zRot];
    segRotMat = pagemtimes(RotMat, TribalSys);
end

end

