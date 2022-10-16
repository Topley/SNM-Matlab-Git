% function AllTrials_Concatenated_Pooled_Coherence_BATCH
% function Run_BATCH_DeltaF_MT
clearvars

%%%% Pick files directory - use wildcard t select specific muscles
rootDir = cd;
masterSheet = fullfile(rootDir,'Master Subject Data Sheet.xls');
% subject = 'Coco01';
%subject = 'Coco03_23661_MG_5_v23decomposed';
[allFilenames] = getCocoFiles(masterSheet, 'muscle', 'TA', 'feedback', 'Hold', 'contraction', 'Coco')
fsamp = 2048;
%  joint2Remove = ~contains(filenames, 'Coco03') ;
%  filenames(joint2Remove) = [];

%%%% User needs to adjust these before running script
groupTAZ ={}; groupSolZ ={}; groupMGZ ={}; groupTASolZ ={}; groupTAMGZ ={}; groupMGSolZ ={};
groupTAF ={}; groupMGF ={}; groupSolF ={}; groupTASolF ={}; groupTAMGF ={}; groupMGSolF ={};

for i = 1:3:length(allFilenames)
    filenames = [allFilenames{i};allFilenames{i+1};allFilenames{i+2}];
    
    subject = filenames(1,1:6)
    
    if contains(subject, 'Coco01')
        usermaxTorque = 52;
    elseif contains(subject, 'Coco02')
        usermaxTorque = 53.5;
    elseif contains(subject, 'Coco03')
        usermaxTorque = 45;
    elseif contains(subject, 'Coco04')
        usermaxTorque = 46;
    end
    
    allFiring = [];
    [allTAFiring, allsmoothTA, allMGFiring, allsmoothMG, allSolFiring, allsmoothSol] = Pooled_Coherence_Concatenate_Trials(fullfile(rootDir, subject, 'decomposed'),filenames, fsamp, usermaxTorque);
        TAFiring = [allTAFiring{1};allTAFiring{2};allTAFiring{3}];
        MGFiring = [allMGFiring{1};allMGFiring{2};allMGFiring{3}];
        SolFiring = [allSolFiring{1};allSolFiring{2};allSolFiring{3}];
%     Loadings ={}; fp1Firing = []; fp2Firing =[]; fp3Firing =[];
%     for j = 1:3
%         TAfiring =allTAFiring{j}';
%         smoothTA = allsmoothTA{j}';
%         TAfiring(:,max(TAfiring)==0) = [];
%         smoothTA(:,max(smoothTA)==0) = [];
%         allFiring = TAfiring;
%         allSmoothFiring = smoothTA;
%         
%         if ~isempty(allSolFiring{j})
%             Solfiring = allSolFiring{j}';
%             smoothSol = allsmoothSol{j}';
%             Solfiring(:,max(Solfiring)==0) = [];
%             smoothSol(:,max(smoothSol)==0) = [];
%             allFiring = [allFiring,Solfiring];
%             allSmoothFiring = [allSmoothFiring,smoothSol];
%         end
%         
%         if ~isempty(allMGFiring{j})
%             MGfiring = allMGFiring{j}';
%             smoothMG = allsmoothMG{j}';
%             MGfiring(:,max(MGfiring)==0) = [];
%             smoothMG(:,max(smoothMG)==0) = [];
%             allFiring = [allFiring,MGfiring];
%             allSmoothFiring = [allSmoothFiring,smoothMG];
%         end
        
%         %[coeff,score,latent,tsquared,explained] = pca(allFiring);
%         try
%             [lambda,psi,T,stats,F] = factoran(allSmoothFiring,3, 'rotate','promax', 'maxit',1500);
%         end
%         
%         Loadings{j} = lambda;
        
        %         try
        %             biplot(lambda(1:size(TAfiring,2),:),'LineWidth',2,'MarkerSize',20, 'Color',[0.1216 0.4667 0.7059],'varlabels',num2str((1:size(TAfiring,2))') )
        %         end
        %         hold on
        %         try
        %             %             biplot(lambda(size(TAfiring,2)+1:size(allFiring,2),:),'LineWidth',2,'MarkerSize',20, 'Color','g','varlabels',num2str((size(TAfiring,2)+1:size(allFiring,2))') )
        %                         %   [lambda,psi,T,stats,F] = factoran(testFiring,3, 'rotate','promax', 'maxit',1500);
        %                         %
        % %                          biplot(lambda(1:size(MGFiring,1),:),'LineWidth',2,'MarkerSize',20, 'Color',[0.1725 0.6275 0.1725],'varlabels',num2str((1:size(MGFiring,1))') )
        %                         % hold on
        %                        % biplot(lambda(size(MGfiring,2)+1:size(allFiring,2),:),'LineWidth',1,'MarkerSize',20, 'Color',[0.8392 0.1529 0.1569],'varlabels',num2str((size(MGfiring,2)+1:size(allFiring,2))') )
        %                          biplot(lambda(size(MGfiring,2)+1:size(allFiring,2),:),'LineWidth',1,'MarkerSize',20, 'Color',[0.1725 0.6275 0.1725],'varlabels',num2str((size(MGfiring,2)+1:size(allFiring,2))') )
        %         end
        %         xlabel('Sol Pool Driven');
        %         ylabel('MG Pool Driven');
        %         zlabel('Coactivation Driven');
        
%         fp1 = Loadings{j}(:,1);
%         fp2 = Loadings{j}(:,2);
%         fp3 = Loadings{j}(:,3);
%         
%         fp1Firing = [fp1Firing; allFiring(:,fp1>.7)'];
%         fp2Firing = [fp2Firing; allFiring(:,fp2>.7)'];
%         fp3Firing = [fp3Firing; allFiring(:,fp3>.7)'];
%     end 
        
        %     set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
        %     set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5
        %     saveas(gcf, fullfile('G:\Coco Subject Data\Group figures',['_', subject,'_PF_AllUnits_FactorAnalysis2.pdf']))
        
%         indx = str2num(subject(end));
%         if ~isempty(fp1Firing)
%             try
%                 [fp1Fig,fp1Z,  fp1F, freqbandFP1] = pooledIntraCoherence('Factor1', fp1Firing,fsamp,1);
%                 %set(TAFig, 'Visible', 'on')
%                 print (fp1Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP1_Concatenated_Coherence.pdf']);
%                 close(fp1Fig);
%                 groupfp1Z{indx} = fp1Z;
%                 groupfp1F{indx} = fp1F;
%             end
%         end
%         
%         if ~isempty(fp2Firing)
%             try
%                 [fp2Fig,fp2Z,  fp2F, freqbandFP2] = pooledIntraCoherence('Factor2', fp2Firing,fsamp,1);
%                 %set(TAFig, 'Visible', 'on')
%                 print (fp2Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP2_Concatenated_Coherence.pdf']);
%                 close(fp2Fig);
%                 groupfp2Z{indx} = fp2Z;
%                 groupfp2F{indx} = fp2F;
%             end
%         end
%         
%         if ~isempty(fp3Firing)
%             try
%                 [fp3Fig, fp3Z,  fp3F, freqbandFP3] = pooledIntraCoherence('Factor 3', fp3Firing,fsamp,1);
%                 %set(TAFig, 'Visible', 'on')
%                 print (fp3Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP3_Concatenated_Coherence.pdf']);
%                 close(fp3Fig);
%                 groupfp3Z{indx} = fp3Z;
%                 groupfp3F{indx} = fp3F;
%             end
%         end
%         
%         
%         try
%             [fp12Fig,fp12Z,  fp12F] = pooledInterCoherence('F1', 'F2', fp1Firing, fp2Firing, fsamp,1);
%             %set(fp12Fig, 'Visible', 'on')
%             print (fp12Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP1_FP2_Pooled_Concatenated_Coherence.pdf']);
%             close(fp12Fig);
%             groupfp12Z{indx} = fp12Z;
%             groupfp12F{indx} = fp12F;
%         end
%         
%         try
%             [fp13Fig,fp13Z,  fp13F] = pooledInterCoherence('F1', 'F3', fp1Firing, fp3Firing, fsamp,1);
%             %set(fp12Fig, 'Visible', 'on')
%             print (fp13Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP1_FP3_Pooled_Concatenated_Coherence.pdf']);
%             close(fp13Fig);
%             groupfp13Z{indx} = fp13Z;
%             groupfp13F{indx} = fp13F;
%         end
%         
%         try
%             [fp23Fig,fp23Z,  fp23F] = pooledInterCoherence('F2', 'F3', fp2Firing, fp3Firing, fsamp,1);
%             %set(fp12Fig, 'Visible', 'on')
%             print (fp23Fig,'-dpdf',['G:\Coco Subject Data\Group figures\', '_', subject,'_FP2_FP3_Pooled_Concatenated_Coherence.pdf']);
%             close(fp23Fig);
%             groupfp23Z{indx} = fp23Z;
%             groupfp23F{indx} = fp23F;
%         end
%     end
    
        if ~isempty(TAFiring)
            try
                [TAFig,TAZ,  TAF, freqbandTA] = Pooled_Intramuscular_Coherence('TA', TAFiring,1,fsamp);
                %set(TAFig, 'Visible', 'on')
                print (TAFig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', '_', subject,'_Coco_TA_Concatenated_Coherence.pdf']);
                close(TAFig);
                groupTAZ{indx} = TAZ;
                groupTAF{indx} = TAF;
            end
        end
    
        if ~isempty(SolFiring)
            try
                [SolFig,SolZ,SolF, freqbandSol] = Pooled_Intramuscular_Coherence('Sol', SolFiring,1,fsamp);
    %             set(SolFig, 'Visible', 'on')
                print (SolFig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', '_', subject, '_Coco_Sol_Concatenated_Coherence.pdf']);
                close(SolFig);
                groupSolZ{indx} = SolZ;
                groupSolF{indx} = SolF;
            end
        end
    
        if ~isempty(MGFiring)
            try
                [MGFig,MGZ,MGF, freqbandMG] = Pooled_Intramuscular_Coherence('MG', MGFiring,1,fsamp);
               % set(MGFig, 'Visible', 'on')
                print (MGFig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', '_', subject, '_Coco_MG_Concatenated_Coherence.pdf']);
                close(MGFig);
                groupMGZ{indx} = MGZ;
                groupMGF{indx} = MGF;
            end
        end
        %         disp(['There is likely an error with a variable in ',filenames(1,1:12)]);
        %         continue
    
        if ~isempty(TAFiring) && ~isempty(SolFiring)
            try
                [pool1Fig, pool1Z,  p1F] = Pooled_Intermuscular_Coherence('TA', 'Sol', TAFiring,SolFiring,1,fsamp);
                print (pool1Fig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', subject, '_Coco_TASol_Concatenated_Coherence.pdf']);
                %set(pool1Fig, 'Visible', 'on');
                close(pool1Fig);
                groupTASolZ{indx} = pool1Z;
                groupTASolF{indx} = pool1F;
            end
        end
    
        if ~isempty(allTAFiring) && ~isempty(allMGFiring)
            try
                [pool2Fig, pool2Z, p2F] = Pooled_Intermuscular_Coherence('TA', 'MG', TAFiring,MGFiring,1,fsamp);
                print (pool2Fig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', subject, '_Coco_TAMG_Concatenated_Coherence.pdf']);
               % set(pool2Fig, 'Visible', 'on');
               close(pool2Fig);
                groupTAMGZ{indx} = pool2Z;
                groupTAMGF{indx} = pool2F;
            end
        end
    
        if ~isempty(allMGFiring) && ~isempty(allSolFiring)
            try
                [pool3Fig,pool3Z,  p3F] = Pooled_Intermuscular_Coherence('MG', 'Sol', MGFiring,SolFiring,1,fsamp);
                print (pool3Fig,'-dpdf',['G:\CoCo Hub\Coco Subject Data\Group figures\', subject, '_Coco_MGSol_Concatenated_Coherence.pdf']);
               % set(pool3Fig, 'Visible', 'on');
               close(pool3Fig);
                groupMGSolZ{indx} = pool3Z;
                groupMGSolF{indx} = pool3F;
            end
        end
    
end

% end