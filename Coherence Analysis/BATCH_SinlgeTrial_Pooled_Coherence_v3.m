function BATCH_SinlgeTrial_Pooled_Coherence_v3
% function Run_BATCH_DeltaF_MT
clearvars 

%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls'); 
% subject = 'Coco01';
%subject = 'Coco03_23661_MG_5_v23decomposed';
[filenames] = getCocoFiles(masterSheet, 'muscle', 'Sol', 'feedback', 'Iramp', 'contraction', 'DF')

 joint2Remove = ~contains(filenames, 'Coco03') ;
 filenames(joint2Remove) = [];

%%%% User needs to adjust these before running script


for i = 1:length(filenames)
     subject = filenames{i}(1:6)
    try
         if contains(filenames(i), 'Coco01')
                usermaxTorque = 52;
            elseif contains(filenames(i), 'Coco02')
                usermaxTorque = 53.5;
            elseif contains(filenames(i), 'Coco03')
                usermaxTorque = 45;
            elseif contains(filenames(i), 'Coco04')
                usermaxTorque = 46;
         end
        Pooled_Coherence_Analysis_v6_MT(fullfile(rootDir, subject, 'decomposed',filenames(i)), 2048,1,500, usermaxTorque)
    catch
        disp(['There is likely an error with a variable in ',filenames(i)]);
        continue
    end
end

end