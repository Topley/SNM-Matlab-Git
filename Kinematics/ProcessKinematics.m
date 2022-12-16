function [JointAngle, JointCenter] = ProcessKinematics(side, staticS, Markers, varargin)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
if nargin > 3
    try
        JointAngle = varargin{1};
        JointCenter = varargin{2};
    catch
        JointAngle = varargin{1};
    end
end
limbs = lower(side);

noRMedial = isfield(Markers, 'RMK');
noLMedial = isfield(Markers, 'LMK');
staticMarker = staticS.staticMarker;

[RotMat, ASISBreadth, HipOrigin, JointCenter] = HipRotationMatrix(Markers);

if contains(limbs, 'right')
    if noRMedial == 0
    [~, ~, nonAnan_RThighR, ~] = nonAnan_RMatrices('RThigh', Markers, JointCenter);
    rotatedRMK3d = pagemtimes(nonAnan_RThighR,'transpose', staticMarker.RMK , 'none');
    rotatedRMK = permute(rotatedRMK3d, [3,1,2]);
    Markers.RMK = rotatedRMK + Markers.RKNE;
    end
        
    JointCenter.RKnee = (Markers.RMK + Markers.RKNE)./2;
    
    [RThigh_R, JointCenter] = segmentRotationMatrix('RThigh', Markers, JointCenter);
    RotMat.RThigh = RThigh_R;
    [RHipX, RHipY, RHipZ] = EulerAngles(RotMat.Pelvis, RotMat.RThigh);
    JointAngle.RHipX = RHipX;
    JointAngle.RHipY = RHipY;
    JointAngle.RHipZ = RHipZ;
    
    [RShank_R, JointCenter] = segmentRotationMatrix('RShank', Markers, JointCenter);
    RotMat.RShank = RShank_R;
    [RKneeX, RKneeY, RKneeZ] = EulerAngles(RotMat.RThigh, RotMat.RShank);
    JointAngle.RKneeX = RKneeX;
    JointAngle.RKneeY = RKneeY;
    JointAngle.RKneeZ = RKneeZ;
    
    [RFoot_R, JointCenter] = segmentRotationMatrix('RFoot', Markers, JointCenter);
    RotMat.RFoot = RFoot_R;
    [RAnkleX, RAnkleY, RAnkleZ] = EulerAngles(RotMat.RShank, RotMat.RFoot);
    JointAngle.RAnkleX = RAnkleX;
    JointAngle.RAnkleY = RAnkleY;
    JointAngle.RAnkleZ = RAnkleZ;
    
else
    if noLMedial == 0
    [~, ~, nonAnan_LThighR, ~] = nonAnan_RMatrices('LThigh', Markers, JointCenter);
    rotatedLMK3d = pagemtimes(nonAnan_LThighR,'transpose', staticMarker.LMK, 'transpose');
    rotatedLMK = permute(rotatedLMK3d, [3,1,2]);
    Markers.LMK = rotatedLMK + Markers.LKNE;
    end 
    JointCenter.LKnee = (Markers.LMK + Markers.LKNE)./2;
    
    [LThigh_R, JointCenter] = segmentRotationMatrix('LThigh', Markers, JointCenter);
    RotMat.LThigh = LThigh_R;
    [LHipX, LHipY, LHipZ] = EulerAngles(RotMat.Pelvis, RotMat.LThigh);
    JointAngle.LHipX = LHipX;
    JointAngle.LHipY = LHipY;
    JointAngle.LHipZ = LHipZ;
    
    [LShank_R, JointCenter] = segmentRotationMatrix('LShank', Markers, JointCenter);
    RotMat.LShank = LShank_R;
    [LKneeX, LKneeY, LKneeZ] = EulerAngles(RotMat.LThigh, RotMat.LShank);
    JointAngle.LKneeX = LKneeX;
    JointAngle.LKneeY = LKneeY;
    JointAngle.LKneeZ = LKneeZ;
    
    [LFoot_R, JointCenter] = segmentRotationMatrix('LFoot', Markers, JointCenter);
    RotMat.LFoot = LFoot_R;
    [LAnkleX, LAnkleY, LAnkleZ] = EulerAngles(RotMat.LShank, RotMat.LFoot);
    JointAngle.LAnkleX = LAnkleX;
    JointAngle.LAnkleY = LAnkleY;
    JointAngle.LAnkleZ = LAnkleZ;
end


end

