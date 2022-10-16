function [upsampled_data] = upsampler(original_data,original_Fs,desired_Fs)

factorFs = desired_Fs/original_Fs;

upsampled_data = zeros(round((size(original_data,1))*factorFs),size(original_data,2));
for i4 = 1:size(original_data,2)
    sig = original_data(:,i4);
    sig_el = 1:numel(sig);
    datapoint = size(upsampled_data,1);

    interp_point = linspace(min(sig_el),max(sig_el),datapoint);
    interp_sig = interp1(sig_el,sig,interp_point);
    upsampled_data(:,i4) = interp_sig;
end

time_original=linspace(0,length(original_data)/original_Fs,length(original_data))';
time_upsampled=linspace(0,length(upsampled_data)/desired_Fs,length(upsampled_data))';

% plot(time_original,original_data(:,3));
% hold on;
% plot(time_upsampled,upsampled_data(:,3));
% legend('original','upsampled')
end 