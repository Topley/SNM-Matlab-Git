function [TABinarySpT, TAFiring, MGBinarySpT, MGFiring, SolBinarySpT, SolFiring] = Pooled_Coherence_Multiple_Trials_Base(filePath, allFilenames, trial)
%close all %-except m1Fig m2Fig poolFig;

%%%%  test file/variables for batch %%%%

%  filename = 'CoCotest02_60009_TA_1p_v23decomposed_MUCLEANED.mat';
%  Muscle1 = 'TA';
%  Muscle2 = 'Sol';
%  Muscle3 = 'MG';
%  N = 10;
%  fsamp = 2048;
%  LW = 1;
%  LimFreq = 500;
%  CoCo = 1;
% coherFolder = 'F:\Toolkit\Mirror\Coherence\DF Trials\Torque Feedback';

%%%%%
warning('off','all')

TAFiring = {}; MGFiring = {}; SolFiring = {};
TABinarySpT = {}; MGBinarySpT = {}; SolBinarySpT = {};
for i = 1:size(allFilenames,1)
    [fileDir,mucleanFile, decompFile,clusterFile, pdfdir] = Set_Directories(fullfile(filePath,allFilenames(i,:)), trial);
    under = strfind(clusterFile,'_');
    decompFileTA = [decompFile(1:under(2)), 'TA',decompFile(under(2):end)];
    clusterFileTA = [clusterFile(1:under(2)), 'TA',clusterFile(under(2):end)];
    decompFileMG = [decompFile(1:under(2)), 'MG',decompFile(under(2):end)];
    clusterFileMG = [clusterFile(1:under(2)), 'MG',clusterFile(under(2):end)];
    decompFileSol = [decompFile(1:under(2)), 'Sol',decompFile(under(2):end)];
    clusterFileSol = [clusterFile(1:under(2)), 'Sol',clusterFile(under(2):end)];
    try
        [TABinarySpT{i}, TAFiring{i}] = Get_Units(fileDir, decompFileTA, clusterFileTA, trial);
    end
    try
        [MGBinarySpT{i}, MGFiring{i}] = Get_Units(fileDir, decompFileMG, clusterFileMG, trial);
    end
    try
        [SolBinarySpT{i}, SolFiring{i}] = Get_Units(fileDir, decompFileSol, clusterFileSol, trial);
    end
end

end

function [binarySpT, binaryFiring] = Get_Units(fileDir, decompFile, clusterFile, trial)
fsamp = 2048;

dVars = matfile(fullfile(fileDir,decompFile));
cVars = matfile(fullfile(fileDir,clusterFile));

MUPulses2 = dVars.MUPulses;
EMGFeedback = filterEMG(cVars.EMG);

%%% Find the contraction hold window based on the EMG activity of muscle
%%% used for visual feedback

if contains(trial, 'hold')
    TraceFeedback = cVars.TraceFeedback;
    [stax] = findHolds(cVars.TraceFeedback,EMGFeedback, 1);
    refSignal = cVars.JR3Mz;
else
    stax = 5;
    refSignal = EMGFeedback;
end

if contains(decompFile, 'Coco01')
    endax = stax+30;
elseif contains(decompFile, 'Coco')
    stax = stax+5;
    endax = stax+60;
else
    stax = stax+5;
    endax = stax+20;
end

try
    [MUST, MUFiring] = binarySpikeTrain(MUPulses2,refSignal );
    [~, binarySpT,binaryFiring, COV, rSpikes] = Remove_MUs_Auto(MUPulses2,MUST,[stax endax]);
end

end

