() = CreatesubjectInformationSheet(subjectfolder)
clear all; clc

rootDir = cd;
subDir = fullfile(rootDir,'Coco04','to_cluster');

filenames = dir(fullfile(subDir, '*.mat'));

Trials = {filenames(:).name};
T = cell2table(Trials')

writetable(T,['Master Subject Data Sheet.xls'],'WriteMode','Append')

           
        