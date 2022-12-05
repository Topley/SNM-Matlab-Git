clear all
close all
clc
%% Delta F Batch run

% Must be in outermost directory containing all study subject data
SubjectID = 'TRDTest04';
SubjectDir = cd;
rootDir = fullfile(SubjectDir, SubjectID);
cleanFiles = dir(fullfile(rootDir,'decomposed', '*MG*MUCLEANED.mat'));
%%%% User needs to adjust these before running script


saveData = 1;              % Save excel file 1 or 0
Auto = [1,0];                   % 1; %no bad MUs 1 or 0                  % if you want to display PDFs 1=summary, 2=summary+units, 3=summary+units+delta f's, 4=all
loop = 0;
usermaxTorque = 50;
groupDF = [];
for i = 1:size(cleanFiles,1)
    try
    filename = fullfile(cleanFiles(i).folder, cleanFiles(i).name);
    [deltaFData] =  delta_f_analysis_TRD_v2(filename, usermaxTorque, Auto, saveData);
    groupDF = [groupDF;deltaFData];
    catch
    end 
end

% end