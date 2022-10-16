%%%% BATCH Process DeltaF's
clearvars % better than clear all

%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
filenames = dir('MartinTest\decomposed\*MUCLEANED.mat')
% filenames(joint2Remove) = [];

%%%% User needs to adjust these before running script


saveexcel =0;              % Save excel file 1 or 0
Auto = 0;                   % 1; %no bad MUs 1 or 0                  % if you want to display PDFs 1=summary, 2=summary+units, 3=summary+units+delta f's, 4=all
keepPics = [1, 2];               % save 1=summary, 2=summary/units, 3=summary/units/deltas
rere = [];
DeRe = [];
peakTorq = [];
loop = 0;
for i = 1:size(filenames,1)

       usermaxTorque = 100;
       ij = 1
       

       [recruitmentThresh,dere] =  delta_f_analysis_TRD_v1(fullfile(filenames(i).folder,filenames(i).name),usermaxTorque,Auto,saveexcel, keepPics, ij);
    
end

% end