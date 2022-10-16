t2 = tiledlayout(1,2)
title(t2,'100 Random CST Permutations vs All Unique Pairs Coherence');
filenames = dir('*coher_v2*.mat');
count = 0;
avg = 0;
pavg = 0;

avgCOF =0;
avgCOF2 = 0;

pavgCOF =0;
pavgCOF2 = 0;

Tableaumap = Tableau;
for i = 1:length(filenames)
  try
    load(filenames(i).name,'rF', 'rz', 'rCOFz', 'rCOF2z', 'pF', 'pz', 'pCOFz', 'pCOF2z')
    nexttile(1)
    hold on
    plot(rF,rz, 'Color', Tableaumap(i,:))
    hold off
    avg = avg+rz;
    avgCOF = avgCOF+rCOFz;
    avgCOF2 = avgCOF2+rCOF2z;
    title('Random Permutations')

   nexttile(2)
    hold on
    plot(pF,pz, 'Color', Tableaumap(i,:))
    hold off
    title('All Unique Pairs')
    pavg = pavg+pz;
    pavgCOF = pavgCOF+pCOFz;
    pavgCOF2 = pavgCOF2+pCOF2z;
    count=count+1;
  catch
  end 
end 

    
    
 avg=avg./count;
 avgCOF = avgCOF./count;
 avgCOF2 = avgCOF2./count;
 
 pavg=pavg./count;
 pavgCOF = pavgCOF./count;
 pavgCOF2 = pavgCOF2./count;

nexttile(1)
hold on
plot(rF,avg,'Color', 'k', 'LineWidth', 2)
plot(rF,ones(1,length(rF))*avgCOF,'r');
plot(rF,ones(1,length(rF))*avgCOF2,'k');
hold off
xlim([0 80])
ylim([0 15])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')

nexttile(2)
hold on
plot(pF,pavg,'Color', 'k', 'LineWidth', 2)
plot(pF,ones(1,length(pF))*pavgCOF,'r');
plot(pF,ones(1,length(pF))*pavgCOF2,'k');
hold off
xlim([0 80])
ylim([0 15])
xlabel('Frequency (Hz)')
ylabel('Coherence (z-transform)')

print (gcf,'-dpdf',['RandPerm_vs_AllPairs_MG_Coherence.pdf']);