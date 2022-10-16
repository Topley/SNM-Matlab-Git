function [poolFig,poolZ,  pF] = pooledInterCoherence(muscle1, muscle2, muscle1Firing,muscle2Firing,LW,fsamp)
% This function will take the smoothed binary spike trains of the active
% "good" units from TWO muscles and run the pooled coherence.
% The autospectra of each unit is calculated from both muscles, and so is the coherence
% spectra of the two separate units within and between muscles.
% The autospectra of each unit is then subtracted out of the pooled
% coherence similar to the single muscle coherence
% wcoherence(sum(muscle1Firing,1), sum(muscle2Firing,1), fsamp)
%%%% summed CST intermuscular coherence
LW = 1;
fsamp = 2048;
% 
%   disp(['Processing ', muscle1, muscle2, ' CST Coherence'])
% try
%     [M1xx,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle1Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%     [M2yy,cstF] = cpsd(detrend(sum(muscle2Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%     [M12xy,cstF] = cpsd(detrend(sum(muscle1Firing,1),0),detrend(sum(muscle2Firing,1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%     
%     cstCOHT = abs(M12xy).^2./(M1xx.*M2yy);
% catch
%     disp('problem with CST coherence')
% end

%%%%%%
% Pooled coherence between agonist and antagonist
%%%%%%
try
    poolFig = figure('Visible','off');
    poolax = axes('Parent',poolFig);
    hold(poolax,'on');
    PxxPool = 0;
    PyyPool = 0;
    PxyPool = 0;
    Unitz = [];
    unitCOHT = [];
    loopCOHT = 0;
    Iterp = 0;
    disp(['Processing ', muscle1, muscle2, ' Shared Residual Coherence'])
    totalSegmentsP = 0;
    for j = 1:size(muscle1Firing,1)
        
        for ij = 1:size(muscle2Firing,1)
            [Pxx,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle1Firing(j,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
            [Pyy,pF] = cpsd(detrend(muscle2Firing(ij,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
            [Pxy,pF] = cpsd(detrend(muscle1Firing(j,:),0),detrend(muscle2Firing(ij,:),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
            PxxPool = PxxPool + Pxx;
            PyyPool = PyyPool + Pyy;
            PxyPool = PxyPool + Pxy;
            unitCOHT = abs(Pxy).^2./((Pxx.*Pyy)+eps);
            loopCOHT = loopCOHT+unitCOHT;
            segments  = floor(size(muscle1Firing,2)/(LW*fsamp));
            totalSegmentsP = totalSegments + segments;
            Unitz = (atanh(sqrt(unitCOHT))) / sqrt(0.5/(segments)); % this is correct
            plot(poolax,pF,Unitz, 'Color', [.7 .7 .7])
            Iterp = Iterp+1;
        end
    end
    
    
    pCOHT = abs(PxyPool).^2./(PxxPool.*PyyPool);
    pCL = 1-0.05^(1/(totalSegmentsP/(fsamp*LW)));
    poolZ = (atanh(sqrt(pCOHT))) / sqrt(0.5/(totalSegments/Iterp)); % this is correct
    biasPool = mean(poolZ(pF>250&pF<500));
%     poolZ = poolZ-biasPool;
%     pCL = pCL-biasPool;
    noisePoolCOFz = max(poolZ(pF>500));
    
  
    plot(poolax,pF,ones(1,length(pF))*pCL,'r', 'LineWidth', 2);
    plot(poolax,pF,ones(1,length(pF))*noisePoolCOFz,'b', 'LineWidth', 2);
    plot(poolax, pF,poolZ, 'Color','k', 'LineWidth', 2)
    xlim(poolax,[0 50]);
    ylim(poolax,[0 20]);
    yticks(poolax,[0 5 10 15 20]);
    xlabel(poolax,'Frequency (Hz)')
    ylabel(poolax,'Coherence (z-transform)')
    title(poolax,[muscle1,'-',muscle2, ' Coherence'])
    hold(poolax,'off');
catch
    disp('problem with pooled coherence')
end
% set(poolFig,'Visible','on');

end