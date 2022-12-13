files = dir('*.mat');

for i = 1:length(files)
    load(files(i).name);
    
    try
        [HipRotMat, ~, ~, JointCenter] = HipRotationMatrix(Markers);
    end
    
    try
        [rTHIRotMat, JointCenter] = segmentRotationMatrix('RThigh', Markers, JointCenter);
    end
    
    try
        [rShinRotMat, JointCenter] = segmentRotationMatrix('RShank', Markers, JointCenter);
    end
    
    try
        [rFootRotMat, JointCenter] = segmentRotationMatrix('RFoot', Markers, JointCenter);
    end
    
    try
        [JointAngle.RHip] = EulerAngles(HipRotMat,rTHIRotMat);
    end
    
    try
        [JointAngle.RKnee] = EulerAngles(rTHIRotMat,rShinRotMat);
    end
    
    try
        [JointAngle.RAnkle] = EulerAngles(rShinRotMat,rFootRotMat);
    end
    
    try
        try
            save(files(i).name, 'JointAngle', 'JointCenter', '-append');
        catch
            
            try
                save(files(i).name, 'JointAngle', '-append');
            catch
                save(files(i).name, 'JointCenter', '-append');
            end
        end
    end
    clearvars -except i files
end
