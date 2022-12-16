function [staticMedialMarker, localLateralMarker, segmentR, JointCenter] = nonAnan_RMatrices(segment, Markers, JointCenter)
%%%% Rotation Matrix
%   First argument specifies the segment the rotation matrix represents.
%   The second argument is a cell array of the markers to build the matrix
%   from
TribalSys = [0,1,0; 1,0,0;0,0,1];
segmentName = segment;
    
switch segmentName
    case 'RThigh'
        planarMarker = Markers.RTHI;
        lateralMarker = Markers.RKNE;
        medialMarker = Markers.RMK;
        proximalMarker = JointCenter.RHip;
        distalMarker = lateralMarker;
        
    case 'LThigh'
        planarMarker = Markers.LTHI;
        lateralMarker = Markers.LKNE;
        medialMarker = Markers.LMK;
        proximalMarker = JointCenter.LHip;
        distalMarker = lateralMarker;
        
    case 'RShank'
        planarMarker = Markers.RTIB;
        lateralMarker = Markers.RANK;
        medialMarker = Markers.RMA;
        proximalMarker = Markers.RKNE;
        distalMarker = lateralMarker;
        
    case 'LShank'
        planarMarker = Markers.LTIB;
        lateralMarker = Markers.LANK;
        medialMarker = Markers.LMA;
        proximalMarker = Markers.LKNE;
        distalMarker = lateralMarker;
        
    otherwise
        disp('No segment chosen')
end

    axisY = proximalMarker - distalMarker;
    tempVect = planarMarker - distalMarker;
    axisX = cross(axisY, tempVect, 2);
    axisZ = cross(axisX, axisY, 2);
    
    Xvect = normr(axisX);
    Yvect = normr(axisY);
    Zvect = normr(axisZ);

    xRot = reshape(Xvect,3, 1, []);
    yRot = reshape(Yvect,3, 1, []);
    zRot = reshape(Zvect,3, 1, []);
    RotMat = [xRot,yRot,zRot];
    segmentR = pagemtimes(RotMat , TribalSys); 
    segmentR = RotMat; 
    
    medialVect3D = medialMarker - distalMarker;
    medial3D = reshape(medialVect3D', [3, 1, size(medialVect3D,1)]);
    rotatedMedialMarker = pagemtimes(segmentR, 'transpose',medial3D, 'none');
    medialVect = reshape(pagetranspose(rotatedMedialMarker),[],3);
    staticMedialMarker = mean(medialVect,1)';
   
    localLateralMarker = lateralMarker - lateralMarker;

end

