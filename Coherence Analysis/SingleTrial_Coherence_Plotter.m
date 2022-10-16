t2 = tiledlayout(1,3)
title(t2,'Group Inter & Intra-muscular Coherence');
% filenames = dir('*coher_v2*.mat');
filenames = dir('*CoCo_coher*.mat');
% filenames = dir('CoCotest02_60010_TA_1p_v23decomposed_MUCLEANED_CoCo_coher_18_31.mat');
count = 0;
avg = 0;
agavg = 0;

avgCOF = 0;
agavgCOF =0;
avgCOF2 = 0;

agavgCOF2 =0;

Tableaumap = Tableau;
for i = 1:length(filenames)
  try
    load(filenames(i).name,'antF', 'antz', 'agF', 'agz', 'agCOFz', 'agCOF2z', 'antCOFz', 'antCOF2z','pF', 'pz', 'pCOFz', 'pCOF2z')
    nexttile(1)
    hold on
    plot(agF,agz, 'Color', Tableaumap(i,:))
    hold off
    title('TA-TA')
    agavg = agavg+agz;
    agavgCOF = agavgCOF+agCOFz;
    agavgCOF2 = agavgCOF2+agCOF2z;

   nexttile(2)
    hold on
    plot(antF,antz, 'Color', Tableaumap(i,:))
    hold off
    title('MG-MG')
    
    nexttile(3)
    hold on
    plot(pF,pz, 'Color', Tableaumap(i,:))
    hold off
    title('TA-MG')
    avg = avg+pz;
    avgCOF = avgCOF+pCOFz;
    avgCOF2 = avgCOF2+pCOF2z;
    count = count+1
  catch
  end 
end 

    avg = avg/count;
    avgCOF = avgCOF/count;
    avgCOF2 = avgCOF2/count;
    agavg = agavg/count;
    agavgCOF = agavgCOF/count;
    agavgCOF2 = agavgCOF2/count;
   
nexttile(1)
hold on
plot(agF,agavg,'Color', 'k', 'LineWidth', 2)
plot(agF,ones(1,length(agF))*agavgCOF,'r');
plot(agF,ones(1,length(agF))*agavgCOF2,'k');
hold off
xlim([0 80])
ylim([0 15])
t2.XLabel.String = 'Frequency (Hz)';
t2.YLabel.String = 'Coherence (z-transform)';

nexttile(2)
hold on
plot(antF,ones(1,length(antF))*antCOFz,'r');
plot(antF,ones(1,length(antF))*antCOF2z,'k');
hold off
xlim([0 80])
ylim([0 15])

nexttile(3)
hold on
plot(pF,avg,'Color', 'k', 'LineWidth', 2)
plot(pF,ones(1,length(pF))*avgCOF,'r');
plot(pF,ones(1,length(pF))*avgCOF2,'k');
hold off
xlim([0 80])
ylim([0 15])


print (gcf,'-dpdf',['TA_and_MG_Intramuscular_Coherence.pdf']);