[num,txt,raw] = xlsread('G:\TRD Hub\Test Subject Data\MartinTest\quietstance_04.csv');
load('MartinTest_25804_Sol_5.mat');

[p,q] = rat(2048/960);

trdForceX = num(5:end,9);

downEMG = num(5:end,17);

testFx = resample(trdForceX, p, q);

upEMG = resample(downEMG, p, q);

[CtrdOTB, LagtrdOTB] = xcorr(CrossCorrEMG,downEMG);

CtrdOTB = CtrdOTB/max(CtrdOTB);

[M21,I21] = max(CtrdOTB);
t21 = LagtrdOTB(I21);

plot(LagtrdOTB,CtrdOTB,[t21 t21],[-0.5 1],'r:')
text(t21+100,0.5,['Lag: ' int2str(t21)])
ylabel('C_{21}')
% axis tight
title('Cross-Correlations')

upEMG2 = downEMG(t21*2:end);
CrossCorrEMG2 = CrossCorrEMG();
figure
ax(1) = subplot(3,1,1);
plot(upEMG2)
ylabel('s_1')
axis tight

ax(2) = subplot(3,1,2);
plot(CrossCorrEMG)
ylabel('s_2')
axis tight

upEMG = resample(upEMG2, p, q);
