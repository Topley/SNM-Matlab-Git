clear all
close all
clc
%% First section needs to be run "twice"
% first with supdata variable uncommented in the modified_suphold dir
% secondly with flexdata variable uncommented in the modified_suphold dir
% "flexdata" & "supdata" are the ONLY variables needed to plot the rest of the script

filenames = dir('*coher*FN.mat');
allexpnames = [];
for i = 1:length(filenames)
    allexpnames{i} = (filenames(i).name(1:8));
end
expnames = unique(allexpnames);
for i = 1:length(expnames)
    loopnames = dir([expnames{i} '*coher*FN.mat']);
    loopz = [];
    for j = 1:length(loopnames)
        load(loopnames(j).name)
        if MUNumber>2            
            loopz = [loopz, ((z))];
        end
    end
       
%     supdata{i,1} = expnames{i};   %these need to be manually (un)commented when folder change from flex to sup
%     supdata{i,2} = mean(loopz,2); %these need to be manually (un)commented when folder change from flex to sup
%     supdata{i,2}(supdata{i,2}<=0) = 0; %these need to be manually (un)commented when folder change from flex to sup

  
    flexdata{i,1} = expnames{i}; %these need to be manually (un)commented when folder change from flex to sup
    flexdata{i,2} = mean(loopz,2);%these need to be manually (un)commented when folder change from flex to sup
    flexdata{i,2}(flexdata{i,2}<=0) = 0;%these need to be manually (un)commented when folder change from flex to sup


end


 clearvars -except F flexdata supdata 
 
 %%
 %Helpful Definitions
%  COHT = coherence function
%  COF = standard confidence interval
%  COF2 = peak COHT > 500
%  z = z transformed COHT
%  COFz = standard z tansform confidence interval
%  COF2z = peak z transform COHT > 500   
c=1
for i=1:length(supdata)% create new variables for flex and sup data that is easier to use
           try
        loopsupdata(c,:) = supdata{i,2}; 
        c = c+1;
    end
end
loopsupdata(loopsupdata<0) = 0;
c=1
for i=1:length(flexdata) 
     try
    loopflexdata(c,:) = flexdata{i,2};
    c = c+1;
    end
end
loopflexdata(loopflexdata<0) = 0;
% counts should = number of subjects
counts = [1,2,3,4,5,6,7,8,9,10,11;...
            1,2,3,4,5,6,7,8,9,10,11];

        tableau =     [0.1216    0.4667    0.7059
    0.6824    0.7804    0.9098
    1.0000    0.4980    0.0549
    1.0000    0.7333    0.4706
    0.1725    0.6275    0.1725
    0.5961    0.8745    0.5412
    0.8392    0.1529    0.1569
    1.0000    0.5961    0.5882
    0.5804    0.4039    0.7412
    0.7725    0.6902    0.8353
    0.5490    0.3373    0.2941
    0.7686    0.6118    0.5804
    0.8902    0.4667    0.7608
    0.9686    0.7137    0.8235
    0.4980    0.4980    0.4980
    0.7804    0.7804    0.7804
    0.7373    0.7412    0.1333
    0.8588    0.8588    0.5529
    0.0902    0.7451    0.8118
    0.6196    0.8549    0.8980];

       
for i = 1:length(counts)
    
    figure(1);hold all
    plot(F,supdata{counts(2,i),2},'Color', tableau(4,:));
    xlim([0 100])
    title('Supination, n=11')
    movegui('northwest');
    allsup(i,:)=supdata{counts(2,i),2};
    
    figure(2);hold all
    plot(F,flexdata{counts(1,i),2}, 'Color', tableau(2,:) );
    xlim([0 100])
    title('Flexion, n=11')
    movegui('northeast');
    allflex(i,:)=flexdata{counts(1,i),2};
    
    figure(3); hold all
    plot(F,flexdata{counts(1,i),2}-supdata{counts(2,i),2})
    xlim([0 100])
    allsub(i,:) = flexdata{counts(1,i),2}-supdata{counts(2,i),2};

end

figure(1)
plot(F,mean(allsup),'Color', tableau(3,:),'LineWidth',2)
xlim([0 100])
figure(2)
plot(F,mean(allflex),'Color', tableau(1,:),'LineWidth',2)
xlim([0 100])
figure(3)
plot(F,mean(allsub),'k','LineWidth',2)
title('AVG DIFFERENCE in Coherence, n=11')
xlabel('Frequency (Hz)') 
ylabel('Coherence (z)') 
xlim([0 100])

 
figure; hold all
plot(F,mean(allsup),'Color', tableau(3,:),'LineWidth',2)
plot(F,mean(allflex),'Color', tableau(1,:),'LineWidth',2)
xlabel('Frequency (Hz)') 
ylabel('Coherence (z)') 
title('AVG Coherence (Flex vs Sup), n=11')
xlim([0 100])

% difference between avg coherence levels (seems to even out around 50 Hz)
figure(9); hold all
bar(mean(allsup(:,F>0 & F<20),2) - mean(allflex(:,F>0 & F<20),2))
title('AVG DIFFERENCE in Coherence btwn 0 - 20 Hz')

[h,p,ci,stats] = ttest(mean(allsup(:,F>0 & F<20),2),mean(allflex(:,F>0 & F<20),2))


 
