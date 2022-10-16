function [true_fsamp, Time_diff] = emg_sync(emg_fast,emg_slow,Fs_fast,Fs_slow,chart)

% inputs:
% emg_fast = common emg signal that was sampled at the higher rate
% emg_slow = common emg signal that was sampled at the lower rate
% Fs_fast = original sampling rate of emg_fast
% Fs_slow = orignal sampling rate of emg_slow
% chart = optional argument for generating a plot of time-aligned data; 0 =
% no plot; 1 = plot

% outputs:
% true_fsamp = the 'true' sampling rate of emg_slow, assuming that emg_fast
% was sampled exactly at its original sampling rate
% Time_diff = how much later emg_slow recording was started relative to
% emg_fast

% % to test
% emg_fast=xcorr_emg_otb;
% emg_slow=xcorr_emg_vicon;
% Fs_fast=2048;
% Fs_slow=960;
% chart=1;

% transposes data if needed
if size(emg_fast,2)>size(emg_fast,1)
    emg_fast=emg_fast';
end
if size(emg_slow,2)>size(emg_slow,1)
    emg_slow=emg_slow';
end

% debiases and normalizes the emg signals
emg_fast=emg_fast-mean(emg_fast);
emg_fast=emg_fast/max(emg_fast);

emg_slow=emg_slow-mean(emg_slow);  
emg_slow=emg_slow/max(emg_slow);

% autocorrelation operation
autocorrelation_matrix = [];

% sets search range for cross-correlation
Fs_min=round(Fs_slow*0.98);
Fs_max=round(Fs_slow*1.02);
new_Fs=Fs_min:0.1:Fs_max;

% cross-correlation operation
for j=new_Fs
    a2 = resample(emg_fast,round(j*5),round(Fs_fast*5)); 
    [acor,lag] = xcorr(emg_slow,a2); %
    autocorrelation_matrix = [autocorrelation_matrix max(abs(acor))];
end
[~,I] = max(autocorrelation_matrix);
true_fsamp = new_Fs(I);

% upsamples slow emg from its true sampling rate (true_fsamp) to the fast
% emg sampling rate and finds the lag between the signals
slow_emg_upsampled = resample(emg_slow,round(Fs_fast*5),round(true_fsamp*5));
[acor,lag] = xcorr(emg_fast,slow_emg_upsampled); %
[~,I] = max(abs(acor));
Time_diff = lag(I)/Fs_fast;


% optional plot
if chart==1
    t_fast=(1:length(emg_fast))/Fs_fast;
    plot(t_fast,emg_fast,'b');hold on;
    t_slow = (1:length(slow_emg_upsampled))/Fs_fast+Time_diff;
    plot(t_slow,slow_emg_upsampled,'r')
    legend('fast','slow')
end