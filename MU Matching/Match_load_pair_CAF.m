%function Match_load_pair_CAF(emgall_1,decomposed_1,emgall_2,decomposed_2,contractionTYPE)
fsamp = 2048;

% emgall_1 = 'Coco01_24931_TA_5.mat';
% decomposed_1 = 'Coco01_24931_TA_5_v23decomposed_MUCLEANED.mat';
% 
% emgall_2 = 'Coco01_24941_TA_5.mat';
% decomposed_2 = 'Coco01_24941_TA_5_v23decomposed_MUCLEANED.mat';

[decomposedFile1,filePath] = uigetfile;
emgall_1 = fullfile(filePath,[decomposedFile1(1:17),'.mat']);
decomposedFile1 = fullfile(filePath,decomposedFile1);

[decomposedFile2, ~] = uigetfile(filePath);
emgall_2 = fullfile(filePath,[decomposedFile2(1:17),'.mat']);
decomposedFile2 = fullfile(filePath,decomposedFile2);



splitname1 = strsplit(emgall_1,'_');
splitname2 = strsplit(emgall_2,'_');

filename = [splitname1{1} '_' splitname1{2} '_' splitname2{2} '_' splitname1{3} '_DF' '_MATCHED_' splitname1{4}];

%filename = 'test_save_name';

pchanNM2(1,:)=[nan 1 2 3 4];
pchanNM2(2,:)=[16 8 7 6 5];
pchanNM2(3,:)=[15 14 13 12 11];
pchanNM2(4,:)=[19 18 17 9 10];
pchanNM2(5,:)=[20 21 22 23 24];
pchanNM2(6,:)=[28 29 30 31 32];
pchanNM2(7,:)=[25 26 27 33 34];
pchanNM2(8,:)=[35 36 37 38 39];
pchanNM2(9,:)=[45 46 47 48 40];
pchanNM2(10,:)=[44 43 42 41 49];
pchanNM2(11,:)=[54 53 52 51 50];
pchanNM2(12,:)=[55 56 64 63 62];
pchanNM2(13,:)=[57 58 59 60 61];

load(emgall_1, 'EMGall');
load(decomposedFile1, 'MUPulses');

%EMGall = diff(EMGall);
[MATRIX1] = EMGreconstruction(EMGall);

for R = 1:12
    for C = 1:5
        SIG{R,C} = squeeze(MATRIX1(R,C,:))';
    end
end    
SIG{1,1} = (SIG{1,2} + SIG{2,2} + SIG{2,1})/3;

% MUPulses(cellfun(@isempty,MUPulses)) = [];
% 
% MUPulses1 = MUPulses;
MUPulses1 = SortUnits(MUPulses);

SIG1 = SIG;

%%%
load(emgall_2, 'EMGall');
load(decomposedFile2, 'MUPulses');

%EMGall = diff(EMGall);
[MATRIX1] = EMGreconstruction(EMGall);


for R = 1:12
    for C = 1:5
        SIG{R,C} = squeeze(MATRIX1(R,C,:))';
    end
end    
SIG{1,1} = (SIG{1,2} + SIG{2,2} + SIG{2,1})/3;

% MUPulses(cellfun(@isempty,MUPulses)) = [];
% 
% MUPulses2 = MUPulses;
MUPulses2 = SortUnits(MUPulses);

SIG2 = SIG;

MATCHED_MU = match_CHRIS_CAF( filename,SIG1,MUPulses1,SIG2,MUPulses2,35e-3,fsamp,'MONO',0.8,10,20,500,60);

% save(filename,'MATCHED_MU')
close all 



