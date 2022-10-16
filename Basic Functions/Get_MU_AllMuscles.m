function [TAFiring,SolFiring, MGFiring, TAEMGFeedback,SolEMGFeedback,MGEMGFeedback] = Get_MU_AllMuscles(filename,clusterFile, coco)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Muscle1 = 'TA';
Muscle2 = 'Sol';
Muscle3 = 'MG';
[fileDir, ~] = fileparts(filename);

% filename = filename{:};
clusterFile = fullfile(fileDir,clusterFile);

if contains(filename, 'TA')
    fileStruct = load(filename, 'MUPulses');
    TAEMGFeedback = load(clusterFile,'EMG');
    TAPulses = fileStruct.MUPulses;
    TAFiring = SortUnits(TAPulses);
    
    m2Filename = replace(clusterFile, Muscle1, Muscle2);
    fileStructM2 = load(m2Filename);
    SolEMGFeedback = filterEMG(fileStructM2.EMG);
    SolFiring = [];
    
    m3Filename = replace(clusterFile, Muscle1, Muscle3);
    fileStructM3 = load(m3Filename);
    MGEMGFeedback = filterEMG(fileStructM3.EMG);
    MGFiring = [];
    
elseif contains(filename, 'Sol')
    fileStruct = load(filename, 'MUPulses');
    SolEMGFeedback = load(clusterFile,'EMG');
    SolPulses = fileStruct.MUPulses;
    SolFiring = SortUnits(SolPulses);
    
    m2Filename = replace(clusterFile, Muscle2, Muscle1);
    fileStructM2 = load(m2Filename);
    TAEMGFeedback = filterEMG(fileStructM2.EMG);
    TAFiring = [];
    
    m3Filename = replace(clusterFile, Muscle2, Muscle3);
    fileStructM3 = load(m3Filename);
    MGEMGFeedback = filterEMG(fileStructM3.EMG);
    MGFiring = [];
    
elseif contains(filename, 'MG')
    fileStruct = load(filename, 'MUPulses');
    MGEMGFeedback = load(clusterFile,'EMG');
    MGPulses = fileStruct.MUPulses;
    MGFiring = SortUnits(MGPulses);
    
    m2Filename = replace(clusterFile, Muscle3, Muscle1);
    fileStructM2 = load(m2Filename);
    TAEMGFeedback = filterEMG(fileStructM2.EMG);
    TAFiring = [];
    
    m3Filename = replace(clusterFile, Muscle3, Muscle2);
    fileStructM3 = load(m3Filename);
    SolEMGFeedback = filterEMG(fileStructM3.EMG);
    SolFiring = [];
end

if contains(filename,'TA') && coco > 0
    m2Filename = replace(filename, Muscle1, Muscle2);
    try
        fileStructM2 = load(m2Filename, 'MUPulses');
        m2Pulses = fileStructM2.MUPulses;
        if length(m2Pulses) > 5
            SolFiring = SortUnits(m2Pulses);
        else
        end
    catch
        disp('No Sol MUCLEANED file for coco');
    end
    
    m3Filename = replace(filename, Muscle1, Muscle3);
    
    try
        fileStructM3 = load(m3Filename, 'MUPulses');
        m3Pulses = fileStructM3.MUPulses;
        if length(m3Pulses) > 5
            MGFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No MG MUCLEANED file');
    end
    
elseif contains(filename,'Sol') && coco > 0
    m2Filename = replace(filename, Muscle2, Muscle1);
    
    try
        fileStructM2 = load(m2Filename, 'MUPulses');
        m2Pulses = fileStructM2.MUPulses;
        if length(m2Pulses) > 5
            TAFiring = SortUnits(m2Pulses);
        else
        end
    catch
        disp('No TA MUCLEANED file for coco');
    end
    
    m3Filename = replace(filename, Muscle2, Muscle3);
    
    try
        fileStructM3 = load(m3Filename, 'MUPulses');
        m3Pulses = fileStructM3.MUPulses;
        if length(m3Pulses) > 5
            MGFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No MG MUCLEANED file');
    end
    
elseif contains(filename,'MG') && coco > 0
    m2Filename = replace(filename, Muscle3, Muscle1);
    
    m3Filename = replace(filename, Muscle2, Muscle3);
    
    try
        fileStructM2 = load(m2Filename, 'MUPulses');
        m2Pulses = fileStructM2.MUPulses;
        if length(m2Pulses) > 5
            TAFiring = SortUnits(m2Pulses);
        else
        end
    catch
    end
    
    try
        fileStructM3 = load(m3Filename, 'MUPulses');
        m3Pulses = fileStructM3.MUPulses;
        if length(m3Pulses) > 5
            SolFiring = SortUnits(m3Pulses);
        else
        end
    catch
        disp('No Sol MUCLEANED file');
    end
    
end

end

