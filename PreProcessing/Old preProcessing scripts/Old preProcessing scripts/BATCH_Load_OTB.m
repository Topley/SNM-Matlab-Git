clear all
close all

% array = {1:64}; %TA

%  array = {65:128}; %MG

% array = {129:192}; %LG

array = {193:256}; %SOL

rootDir = ('G:\TRD Hub\Test Subject Data\TRDTest_02b')
% rootDir = ('G:\CoCo Hub\Coco Subject Data\Coco05')
filenames = dir(fullfile(rootDir,'*.otb'));
c3dFiles = dir(fullfile(rootDir, '*.c3d'));

removeStatic = contains({c3dFiles.name},{'static','Static'});
c3dFiles(removeStatic,:) = [];

stickers = 4;
auxChans = 2;
totalChans = stickers * 64 + auxChans;
    for i = 1:length(filenames)
            fullFilePathOTB = fullfile(filenames(i).folder,filenames(i).name);
            fullFilePathC3D = fullfile(c3dFiles(i).folder,c3dFiles(i).name);
            try
                Load_TRDdata(fullFilePathOTB, fullFilePathC3D, array, totalChans, {''}, 5, 2048)
                %Load_OTBdata_Basic(fullFilePathOTB,array,269,{''},5,2048)
            catch
                disp([filenames(i).name, ' Failed'])
            end
    end
