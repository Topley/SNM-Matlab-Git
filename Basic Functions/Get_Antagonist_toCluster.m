function [AntagEMGFeedback] = Get_Antagonist_toCluster(clusterFile, agonist, antagonist)

antagonistClusFile = [];
fsamp = 2048;

try
    antagonistClusFile = replace(clusterFile, agonist, antagonist);

    load(antagonistClusFile);
    
    goodantEMG = EMG;
    
    AntagEMG = abs(goodantEMG(27,:));
    
    AntagEMG = AntagEMG-mean(AntagEMG(1:fsamp));
    
    antFeedback = rms(AntagEMG, 500);
    
    AntagEMGFeedback = movmedian(antFeedback, fsamp);
    
catch
    disp('Cannot find antagonist cluster file');
end


end

