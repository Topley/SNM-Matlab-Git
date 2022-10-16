function [TAFiring,SolFiring, MGFiring, TAEMGFeedback,SolEMGFeedback,MGEMGFeedback] = getAllCoherenceUnits(decompFile,clusterFile, coco)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Muscle1 = 'TA';
Muscle2 = 'Sol';
Muscle3 = 'MG';
[fileDir, ~] = fileparts(decompFile);

% filename = filename{:};
clusterFile = fullfile(fileDir,clusterFile);

if contains(decompFile, 'TA')
    dVars = matfile(decompFile);
    cVars = matfile(clusterFile);
    TAEMGFeedback = cVars.EMG(27,:);
    TAPulses = dVars.MUPulses;
    TAFiring = SortUnits(TAPulses);
    
    m2Filename = replace(clusterFile, Muscle1, Muscle2);
    cVarsM2 = matfile(m2Filename);
    SolEMGFeedback = filterEMG(cVarsM2.EMG);
    SolFiring = [];
    
    m3Filename = replace(clusterFile, Muscle1, Muscle3);
    cVarsM3 = matfile(m3Filename);
    MGEMGFeedback = filterEMG(cVarsM3.EMG);
    MGFiring = [];
    
elseif contains(decompFile, 'Sol')
    dVars = matfile(decompFile);
    cVars = matfile(clusterFile);
    SolEMGFeedback = clusterFile.EMG(27,:);
    SolPulses = dVars.MUPulses;
    SolFiring = SortUnits(SolPulses);
    
    m2Filename = replace(clusterFile, Muscle2, Muscle1);
    cVarsM2 = matfile(m2Filename);
    TAEMGFeedback = filterEMG(cVarsM2.EMG);
    TAFiring = [];
    
    m3Filename = replace(clusterFile, Muscle2, Muscle3);
    cVarsM3 = matfile(m3Filename);
    MGEMGFeedback = filterEMG(cVarsM3.EMG);
    MGFiring = [];
    
elseif contains(decompFile, 'MG')
    dVars = matfile(decompFile);
    cVars = matfile(clusterFile);
    MGEMGFeedback = cVars.EMG(27,:);
    MGPulses = dVars.MUPulses;
    MGFiring = SortUnits(MGPulses);
    
    m2Filename = replace(clusterFile, Muscle3, Muscle1);
    cVarsM2 = loadmatfile(m2Filename);
    TAEMGFeedback = filterEMG(cVarsM2.EMG);
    TAFiring = [];
    
    m3Filename = replace(clusterFile, Muscle3, Muscle2);
    cVarsM3 = matfile(m3Filename);
    SolEMGFeedback = filterEMG(cVarsM3.EMG);
    SolFiring = [];
end

if contains(decompFile,'TA') && coco > 0
    m2Filename = replace(decompFile, Muscle1, Muscle2);
    try
        dVarsM2 = matfile(m2Filename);
        m2Pulses = dVarsM2.MUPulses;
        if length(m2Pulses) > 5
            SolFiring = SortUnits(m2Pulses);
        else
        end
    catch
        disp('No Sol MUCLEANED file for coco');
    end
    
    m3Filename = replace(decompFile, Muscle1, Muscle3);
    
    try
        dVarsM3 = matfile(m3Filename);
        m3Pulses = dVarsM3.MUPulses;
        if length(m3Pulses) > 5
            MGFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No MG MUCLEANED file');
    end
    
elseif contains(decompFile,'Sol') && coco > 0
    m2Filename = replace(decompFile, Muscle2, Muscle1);
    
    try
        dVarsM2 = matfile(m2Filename);
        m2Pulses = dVarsM2.MUPulses;
        if length(m2Pulses) > 5
            TAFiring = SortUnits(m2Pulses);
        else
        end
    catch
        disp('No TA MUCLEANED file for coco');
    end
    
    m3Filename = replace(decompFile, Muscle2, Muscle3);
    
    try
        dVarsM3 = matfile(m3Filename);
        m3Pulses = dVarsM3.MUPulses;
        if length(m3Pulses) > 5
            MGFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No MG MUCLEANED file');
    end
    
elseif contains(decompFile,'MG') && coco > 0
    m2Filename = replace(decompFile, Muscle3, Muscle1);
    
    m3Filename = replace(decompFile, Muscle2, Muscle3);
    
    try
        dVarsM2 = matfile(m2Filename);
        m2Pulses = dVarsM2.MUPulses;
        if length(m2Pulses) > 5
            TAFiring = SortUnits(m2Pulses);
        else
        end
    catch
    end
    
    try
        dVarsM3 = matfile(m3Filename);
        m3Pulses = dVarsM3.MUPulses;
        if length(m3Pulses) > 5
            SolFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No Sol MUCLEANED file');
    end
    
end

end

