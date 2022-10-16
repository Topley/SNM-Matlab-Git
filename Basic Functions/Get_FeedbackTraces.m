function [clusterFile, Torque, TorqueFeedback, EMGFeedback] = Get_FeedbackTraces(DecomposedFile)
fsamp = 2048;

try
    Under = strfind(DecomposedFile,'_');                      % search for file delimiters
    clusterFile = [DecomposedFile(1:(Under(end-1)-1)),'.mat'];   % Find cluster file if neccessary
catch
    disp('Cannot find Analogs');
end

try
    load(clusterFile);
    
    if exist('AuxChan') == 1
        
        auxTorque = AuxChan(1,:);
        Torque = -(auxTorque-(mean(auxTorque(1:fsamp))));
        
        auxTorqueFeedback = AuxChan(9,:);
        TorqueFeedback = -(auxTorqueFeedback-(mean(auxTorqueFeedback(1:fsamp))));
        EMG = sum(abs(EMG));
        
    else
        
        Torque = (Torque-(mean(Torque(1:fsamp))));
        TorqueFeedback = (JR3Mz-(mean(JR3Mz(1:fsamp))));
        EMG = EMGTorque;
    end
    
catch
    try
        load(DecomposedFile);
        Torque = (Torque-(mean(Torque(1:fsamp))));
        TorqueFeedback = (JR3Mz-(mean(JR3Mz(1:fsamp))));
        EMG = EMGTorque;
    catch
        disp('Cannot find Analogs');
    end
end

try
    if mean(TorqueFeedback)<0
        TorqueFeedback = -TorqueFeedback;
        Torque = -Torque;
    end
    
    absRawEMG = abs(EMG);
    rectifiedEMG = absRawEMG-mean(absRawEMG(1:fsamp));
    rmsEMG = rms(rectifiedEMG, 500);
    EMGFeedback = movmedian(rmsEMG, fsamp);
end

end
