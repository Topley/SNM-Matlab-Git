function [cohFig, cohZ, cohF, cohCL, noiseCLz, freqBands] = Pooled_Intramuscular_Coherence(muscle1, firing, LW, fsamp)

% This function will take the smoothed binary spike trains of the active
% "good" units from ONE muscle and run the pooled coherence.
% The autospectra of each unit is calculated, and so is the coherence
% spectra of the two separate units. The autospectra of each unit is then
% subtracted out of the pooled coherence

% LW = 1;
% fsamp = 2048;

%% Random permutations of active units
% try
%     N = 10;
%     NREAL = size(firing,1);
%     group1 = [1:2:NREAL]; % take odds
%     group2 = [2:2:NREAL]; % take evens
%     if N<round(size(firing,1)/2-1)
%         NI = N;
%     else
%         NI = round(size(firing,1)/2-1);
%     end
%
%     disp(['Processing ', muscle1, ' CST Coherence'])
%     Pxx_CST = 0;
%     Pyy_CST = 0;
%     Pxy_CST = 0;
%     Iter = 100;
%     for t = 1:Iter
%         [Pxx,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group1(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%         [Pyy,cstF] = cpsd(detrend(sum(firing(group2(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%         [Pxy,cstF] = cpsd(detrend(sum(firing(group1(1:NI),:),1),0),detrend(sum(firing(group2(1:NI),:),1),0),hanning(round(LW*fsamp)),0,10*fsamp,fsamp);
%         Pxx_CST = Pxx_CST + Pxx;
%         Pyy_CST = Pyy_CST + Pyy;
%         Pxy_CST = Pxy_CST + Pxy;
%     end
%     cstCOHT = abs(Pxy_CST).^2./(Pxx_CST.*Pyy_CST);
% catch
%     disp('Error processing single muscle CST intra-muscular coherence')
% end

%% All pairs of agonist coherence
try
    cohFig = figure('Visible', 'off');
    cohFax = axes('Parent', cohFig);
    hold(cohFax, 'on');
    pool_Pxx = 0;
    pool_Pyy = 0;
    pool_Pxy = 0;
    unit_Z = [];
    unit_Coh = [];
    loop_Coh = 0;
    Iter = 0;
    disp(['Processing ', muscle1, muscle1, ' Shared Residual Coherence'])
    totalSegments = 0;
    for j = 1:size(firing, 1) - 1
        
        for ij = j + 1:size(firing, 1)
            [MU_Pxx, cohF] = cpsd(detrend(firing(j,:), 0), detrend(firing(j,:), 0), hanning(round(LW * fsamp)), 0, 10 * fsamp, fsamp);
            [MU_Pyy, cohF] = cpsd(detrend(firing(ij,:) ,0), detrend(firing(ij,:), 0), hanning(round(LW * fsamp)), 0, 10 * fsamp, fsamp);
            [MU_Pxy, cohF] = cpsd(detrend(firing(j,:) ,0), detrend(firing(ij,:), 0), hanning(round(LW * fsamp)), 0, 10 * fsamp, fsamp);
            pool_Pxx = pool_Pxx + MU_Pxx;
            pool_Pyy = pool_Pyy + MU_Pyy;
            pool_Pxy = pool_Pxy + MU_Pxy;
            
            unit_Coh = abs(MU_Pxy) .^ 2 ./ (MU_Pxx .* MU_Pyy);
            loop_Coh = loop_Coh + unit_Coh;
            
            unitSegments  = floor(size(firing,2) / (LW * fsamp));
            totalSegments = totalSegments + unitSegments;
            
            unit_Z = (atanh(sqrt(unit_Coh))) / sqrt( 0.5 / (unitSegments)); % this is correct
            plot(cohFax, cohF, unit_Z, 'Color', [.7 .7 .7])
            
            hold on
            Iter = Iter + 1;
        end
    end
    
    poolCoh = abs(pool_Pxy) .^2 ./ (pool_Pxx .* pool_Pyy);
    cohCL = 1 - 0.05 ^ (1 / (totalSegments - 1));
    cohZ = (atanh(sqrt(poolCoh))) / sqrt( 0.5 / (totalSegments / Iter)); % this is correct
    
    freqBands.Delta = trapz(cohF(1:50), cohZ(1:50));    % 0 - 5Hz
    freqBands.Alpha = trapz(cohF(50:150), cohZ(50:150));    % 5 - 15Hz
    freqBands.Beta = trapz(cohF(150:350), cohZ(150:350));    % 15 - 35Hz
    area(cohFax, cohF(1:50), cohZ(1:50), 'SeriesIndex', 1, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(cohFax, cohF(50:150), cohZ(50:150), 'SeriesIndex', 3, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    area(cohFax, cohF(150:350), cohZ(150:350), 'SeriesIndex', 5, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    
    noiseCLz = max(cohZ(cohF>500));
    
    plot(cohFax, cohF, ones(1, length(cohF)) * cohCL, 'r', 'LineWidth', 2);
    plot(cohFax, cohF, ones(1, length(cohF)) * noiseCLz, 'b', 'LineWidth', 2);
    plot(cohFax, cohF, cohZ, 'Color', 'k' , 'LineWidth', 2)
    xlim(cohFax, [0 50]);
    ylim(cohFax, [0 20]);
    yticks(cohFax, [0 5 10]);
    xlabel(cohFax, 'Frequency (Hz)')
    ylabel(cohFax, 'Coherence (z-transform)')
    title(cohFax, [muscle1, '-', muscle1, ' Pooled Coherence'])
    hold(cohFax, 'off');
%     set(cohFig,'Visible','on');
catch
    disp('Error processing single muscle all unique pairs intra-muscular coherence')
end

end