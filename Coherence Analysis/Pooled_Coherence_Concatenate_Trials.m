function [allTAFiring, allsmoothTA, allMGFiring, allsmoothMG, allSolFiring, allsmoothSol] = Pooled_Coherence_Concatenate_Trials(filePath, fullFilenames, fsamp, usermaxTorque)
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
hold1 = fullFilenames(1,:);
hold2 = fullFilenames(2,:);
hold3 = fullFilenames(3,:);

[fileDir,Cleanedfile1,Decompfile1,clusterFile1, pdfdir] = Set_Directories(fullfile(filePath,hold1), 'hold');
[~,Cleanedfile2,Decompfile2,clusterFile2, ~] = Set_Directories(fullfile(filePath,hold2), 'hold');
[~,Cleanedfile3,Decompfile3,clusterFile3, ~] = Set_Directories(fullfile(filePath,hold3), 'hold');

[TAFiring1,smoothTA1, MGFiring1, smoothMG1, SolFiring1, smoothSol1] = oneFileUnits(fileDir, Cleanedfile1, clusterFile1, usermaxTorque);
[TAFiring2, smoothTA2, MGFiring2, smoothMG2, SolFiring2, smoothSol2] = oneFileUnits(fileDir, Cleanedfile2, clusterFile2, usermaxTorque);
[TAFiring3, smoothTA3, MGFiring3, smoothMG3, SolFiring3, smoothSol3] = oneFileUnits(fileDir, Cleanedfile3, clusterFile3, usermaxTorque);

try
allTAFiring{1} = TAFiring1;
allTAFiring{2} = TAFiring2;
allTAFiring{3} = TAFiring3;

allsmoothTA{1} = smoothTA1;
allsmoothTA{2} = smoothTA2;
allsmoothTA{3} = smoothTA3;
end

try
allMGFiring{1} = MGFiring1;%MGFiring2;MGFiring3];
allMGFiring{2} = MGFiring2;
allMGFiring{3} = MGFiring3;

allsmoothMG{1} = smoothMG1;
allsmoothMG{2} = smoothMG2;
allsmoothMG{3} = smoothMG3;
end

try
allSolFiring{1} = SolFiring1;%;SolFiring2;SolFiring3];
allSolFiring{2} = SolFiring2;
allSolFiring{3} = SolFiring3;

allsmoothSol{1} = smoothSol1;
allsmoothSol{2} = smoothSol2;
allsmoothSol{3} = smoothSol3;
end

end

function [TAFiring, smoothTA, MGFiring,smoothMG, SolFiring,smoothSol] = oneFileUnits(fileDir, filename, clusterFile, usermaxTorque);
fsamp = 2048;

varsClean = matfile(fullfile(fileDir,filename));
varsCluster = matfile(fullfile(fileDir,clusterFile));
displayPulses = varsClean.MUPulses;
TraceFeedback = varsCluster.TraceFeedback;
EMGFeedback = EMGpro(varsCluster.EMG, 'filter', {'filtRMS', 500 200});

%%% Find the contraction hold window based on the EMG activity of muscle
%%% used for visual feedback

coco = 1;

[stax] = Set_HoldWindow(TraceFeedback,EMGFeedback);
timeTQ = [0:length(TraceFeedback)-1]/fsamp;

if contains(filename, 'Coco01')
    endax = stax+30;
else
    stax = stax+5;
    endax = stax+60;
end
TAFiring = [];
MGFiring = [];
SolFiring = [];
% If a multiple muscles are being analyzed, load other muscles file and MUs
% select contraction window based on range where motor units from all files
% are active
smoothTA = [];
smoothMG= [];
smoothSol =[];
[TAUnits,SolUnits, MGUnits,TAEMGFeedback,SolEMGFeedback,MGEMGFeedback] = getAllCoherenceUnits(fullfile(fileDir,filename),clusterFile, coco);

if size(TAUnits,2) >4
    try
        [TAST, firingTA] = BinarySpikeTrain(TAUnits, varsCluster.JR3Mz);
        [TAgoodUnits, TAFiring,smoothTA, TACOV, TArSpike] = Remove_MUs_Auto(TAUnits,TAST,[stax endax]);
    end
end

if size(MGUnits,2) >4
    try
        [MGST] = BinarySpikeTrain(MGUnits, varsCluster.JR3Mz);
        [MGgoodUnits,MGFiring, smoothMG,MGCOV, MGrSpike] = Remove_MUs_Auto(MGUnits,MGST,[stax endax]);
    catch
        MGFiring = [];
    end
end

if size(SolUnits,2) >4
    try
        [SolST] = BinarySpikeTrain(SolUnits, varsCluster.JR3Mz);
        [SolgoodUnits,SolFiring,smoothSol, SolCOV, SolrSpike] = Remove_MUs_Auto(SolUnits,SolST,[stax endax]);
    catch
        SolFiring = [];
    end
end

end

