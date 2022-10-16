function [] = Plot_MU_Firing(MUTime,windowST, windowDR,timeTorque,newTorque,newTorqueFeedback, AVERAGE)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;

f3 = figure(3); hold all
xlim([MUTime(1) MUTime(2)]);
ylim([0 25]);
c = 1;
Tableaumap = Tableau;
for i = 1:length(windowDR)
    try
        plot(windowST{i},windowDR{i},'LineWidth',1.75,'SeriesIndex',c);
        c=c+1;
    catch
        c=1;
        plot(windowST{i},windowDR{i},'LineWidth',1.75,'SeriesIndex',c);
    end
    %title({[filename(1:end-4),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],['All ',num2str(length(MUFiring)),' units']},'interpreter','none')
end

plot(timeTorque,newTorque,'k','LineWidth',1.75)
plot(timeTorque,newTorqueFeedback,'k','LineWidth',1.75)
plot([1/fsamp:1/fsamp:size(AVERAGE,2)/fsamp],AVERAGE(size(AVERAGE,1),:),'Color', [.1 .1 .1],'LineWidth',1.75)
%%%
windowDR{size(AVERAGE,1)} = AVERAGE(size(AVERAGE,1),:);
windowST{size(AVERAGE,1)} = 1/fsamp:1/fsamp:(size(AVERAGE,2))/fsamp;

end

