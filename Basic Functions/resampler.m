function [resampled_data] = resampler(original_data,original_Fs,desired_Fs,chart)


% % % to test
% original_data=vicon_data;
% original_Fs=true_fsamp;
% desired_Fs=512;

factorFs = desired_Fs/original_Fs;

resampled_data = zeros(round((size(original_data,1))*factorFs),size(original_data,2));
for i4 = 1:size(original_data,2)
    sig = original_data(:,i4);
    sig_el = 1:numel(sig);
    datapoint = size(resampled_data,1);

    interp_point = linspace(min(sig_el),max(sig_el),datapoint);
    interp_sig = interp1(sig_el,sig,interp_point);
    resampled_data(:,i4) = interp_sig;
end

if chart==1
    time_original=linspace(0,length(original_data)/original_Fs,length(original_data))';
    time_resampled=linspace(0,length(resampled_data)/desired_Fs,length(resampled_data))';
    figure(99)
    plot(time_original,original_data(:,1));
    hold on;
    plot(time_resampled,resampled_data(:,1));
    legend('original','resampled')
    hold off
end