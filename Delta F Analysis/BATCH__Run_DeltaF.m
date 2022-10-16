% function Run_BATCH_DeltaF_MT
clearvars % better than clear all

%  [num,txt,raw] = xlsread('Master Subject Data Sheet.xls') ;
% 
% Subject = 'Coco04';
% filenames = cellfun(@num2str,raw,'un',0);
% % 
% filenames = dir(fullfile(rootDir, '**', '*TA_*.xls'));
%  
% contraction2Remove = contains({filenames.folder}, 'keep') ;
% filenames(contraction2Remove) = [];
% 
% joint2Remove = ~contains({filenames.folder}, 'Torque Ramps') ;
% filenames(joint2Remove) = [];

%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls'); 
% subject = 'Coco01';
%subject = 'Coco03_23661_MG_5_v23decomposed';
[filenames] = getCocoFiles(masterSheet, 'remove', 1,'muscle', 'MG', 'feedback', 'Iramp', 'contraction', 'PF')
% 
 joint2Remove = ~contains(filenames, 'Coco03_23653_MG') ;
 filenames(joint2Remove) = [];

%%%% User needs to adjust these before running script


saveexcel =1;              % Save excel file 1 or 0
Auto = 0;                   % 1; %no bad MUs 1 or 0                  % if you want to display PDFs 1=summary, 2=summary+units, 3=summary+units+delta f's, 4=all
keepPics = [1, 2];               % save 1=summary, 2=summary/units, 3=summary/units/deltas
rere = [];
DeRe = [];
peakTorq = [];
loop = 0;
for i = 1:size(filenames,1)
 subject = filenames{i}(1:6)
    for ij = 1:3
loop = loop+1;
        warning('off', 'MATLAB:LOAD:VariableNotFound'); % disable warning if variable can't be loaded
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
%             fileparts(filenames(i).folder)
       [recruitmentThresh,dere] =  delta_f_CoCoFunct_wSubFunctions_v1(fullfile(rootDir,subject,'decomposed',filenames(i)),usermaxTorque,Auto,saveexcel, keepPics, ij);
          %[recruitmentThresh,dere] =  delta_f_CoCoFunct_v2(fullfile(rootDir,subject,'decomposed',filenames(i)),usermaxTorque,Auto,Pics,saveexcel, keepPics, ij);
%            [recruitmentThresh,dere] =  delta_f_CoCoFunct_v2(fullfile(rootDir,subject,'decomposed',filenames(i)),usermaxTorque,Auto,Pics,saveexcel, keepPics, ij);
        rere= [rere,recruitmentThresh];
        DeRe= [DeRe,dere];
        catch
            disp(['Could not complete delta f for', ' Ramp ', num2str(ij), ' in ', filenames{i}]);
            continue
        end
    end
    
end

% end