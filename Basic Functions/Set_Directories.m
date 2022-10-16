function [subDirectory, trialname, mucleanFile, decompFile, clusterFile, kkFile, pdfdir] = Set_Directories(FullFileName, trialType)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[subDirectory,inputFile] = fileparts(FullFileName);
[subDirectory,~] = fileparts(subDirectory);
under = strfind(inputFile, '_');
trialname = inputFile(1:under(3)-1);

if contains(inputFile, '.mat')
    inputFile = erase(inputFile,'.mat');
end

if  ~contains(inputFile, '_v23decomposed') && length(under) < 3
    clusterFile = fullfile(subDirectory, 'to_cluster', [inputFile, '_5.mat']);
    mucleanFile = fullfile(subDirectory, 'decomposed', [inputFile, '_5_v23decomposed_MUCLEANED.mat']);
    decompFile = fullfile(subDirectory, 'decomposed', [inputFile,'_5_v23decomposed.mat']);
    kkFile = fullfile(subDirectory, 'kk files', [inputFile(1:under(3)),'kk.mat']);
    
elseif  ~contains(inputFile, '_v23decomposed') && length(under) > 3
    clusterFile = fullfile(subDirectory, 'to_cluster', [inputFile(1:under(4)-1), '.mat']);
    mucleanFile = fullfile(subDirectory, 'decomposed', [inputFile, '_v23decomposed_MUCLEANED.mat']);
    decompFile =  fullfile(subDirectory, 'decomposed', [inputFile,'_v23decomposed.mat']);
    kkFile = fullfile(subDirectory, 'kk files', [inputFile(1:under(3)),'kk.mat']);
    
elseif ~contains(inputFile, 'MUCLEANED')
    clusterFile = fullfile(subDirectory, 'to_cluster', [inputFile(1:under(4)-1),'.mat']);
    decompFile =  fullfile(subDirectory, 'decomposed', [inputFile, '_v23decomposed.mat']);
    mucleanFile = fullfile(subDirectory, 'decomposed', [inputFile,'_MUCLEANED.mat']);
    kkFile = fullfile(subDirectory, 'kk files', [inputFile(1:under(3)),'kk.mat']);
else
    clusterFile = [inputFile(1:under(5)-1),'.mat'];
    decompFile = [inputFile(1:under(6)-1),'.mat'];
    mucleanFile = [inputFile,'.mat'];
    kkFile = fullfile(subDirectory, 'kk files', [inputFile(1:under(3)),'kk.mat']);
    
end


switch trialType
    case 'ramp'
        pdfdir = [subDirectory '\deltaf_pdf'];
        if exist(pdfdir,'dir') == 0;
            mkdir([subDirectory '\deltaf_pdf']);
        end
        pdfdir = [subDirectory '\deltaf_pdf'];
    case 'hold'
        pdfdir = [subDirectory '\coherence_pdf'];
        if exist(pdfdir,'dir') == 0;
            mkdir([subDirectory '\coherence_pdf']);
        end
        pdfdir = [subDirectory '\coherence_pdf'];
    case 'TRD'
        pdfdir =[];
end


end

