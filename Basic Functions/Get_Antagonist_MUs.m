function [AntagonistMUs, binaryAntagonistFiring] = Get_Antagonist_MUs(decomposedFile, agonist, antagonist)

antagonistDecomposedFile = [];

try
    antagonistDecomposedFile = replace(decomposedFile, agonist, antagonist);
    load(antagonistDecomposedFile);
    
catch
    disp('Cannot find antagonist decompsoed file');
end

try
    MotorUnits = MUPulses;
    
    [~,FiringLen] = sort(cellfun(@length,MotorUnits));
    MUFiring = MotorUnits(FiringLen);
    
    for MU = 1:size(MUFiring,2)
        index = MUFiring{MU};
        binaryFiring(MU,index) = 1;
    end
    
    AntagonistMUs = MUFiring;
    binaryAntagonistFiring(:,length(binaryFiring):length(trialLength))=0;
    
catch
    AntagonistMUs = [];
    binaryAntagonistFiring=[];
end

end

