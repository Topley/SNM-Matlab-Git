function [rawEMG, filtEMG] = EMGpro(arrayEMG, varargin)
%% EMGpro.v1 - please uplaod any changes to the lab drive with comments
% This function will process single or multiple channels of EMG data from
% the OTB amplifier

%   input argument 1 requires raw EMG data in vector or array format
%   input argument 2 is variable length
%   1st you must declare the channel if you are processing single channel
%   data, if calculating the average across all channels then ignore
%   Any specified processing must be input as a cell array starting with
%   'filter' by the filter type and its parameters.
%   No additional inputs to the EMG variable will output the ABS

%   ex. EMGpro('channel', 27, 'butter', {4, 20/1024, 'high'})
%   20/1024 because 0.5 * fsamp = 1024 with a 20Hz cut
%   creates 4th order 20Hz cutoff high-pass butterworth filter


%% Input argument parsing loop
% variable argument inputs must have a string - the type of processing
% requested - and then based on the processing filter and window properties
% are also required

fsamp = 2048; % always declare this at the beginning of our functions to be safe
intWindow = fsamp/2;
rmsWindow = fsamp/2;
order = 2;
fs = (fsamp/2) / 20;
band = 'low';
if size(varargin,2) > 0
    for i=1:2:3
        switch varargin{i} % parsing inputs
            case 'channel'
                channel = varargin{i+1};
            case 'filter'
                filtType = varargin{i+1};
                try
                    switch filtType
                        case 'rawRMS'
                            rmsL = fsamp/2; %unfiltered RMS length
                        case 'iEMG'
                            intWindow = varargin{i+2};
                        case 'filtRMS'
                            rmsL = fsamp/2; %RMS length
                            rmsWindow = varargin{i+2};
                        case 'butter'
                            order = varargin{i+2}; % filter order
                            fs = (fsamp/2) / varargin{i+3}; % frequency cutoff/band
                            band = varargin{i+4}; % filter type - high, low, bandpass
                    end
                end
            otherwise
                
        end
    end
end

if exist('channel', 'var')
    rawEMG = abs(arrayEMG(channel, :));
else
    rawEMG = abs(sum(arrayEMG)) ./size(arrayEMG, 1);
end

if nargout > 1
    try
        switch filtType
            case 'rawRMS'
                rmsEMG = sqrt(rms(rawEMG .^2, rmsL));
                filtEMG = rmsEMG - mean(rmsEMG(1:fsamp));
            case 'iEMG'
                iEMG = zeros(size(rawEMG));
                for n=1:length(rawEMG)-intWindow
                    if (n>intWindow)
                        iEMG(n) = trapz(rawEMG(n:n+intWindow))/intWindow;
                    end
                end
                frontBuffer = iEMG(intWindow + 1 : intWindow * 2);
                endBuffer = iEMG(end - intWindow * 2 : end - intWindow);
                iEMG(1:length(frontBuffer)) = frontBuffer;
                iEMG(end - length(endBuffer) + 1 : end) = endBuffer;
                filtEMG = iEMG - mean(iEMG(1:fsamp));
            case 'filtRMS'
                filtRMSEMG = movmean(sqrt(rms(rawEMG .^2, rmsL)), rmsWindow);
                filtEMG = filtRMSEMG - mean(filtRMSEMG(1:fsamp));
            case 'butter'
                EMGoutput = sqrt(rms(rawEMG .^2, 500));
                [b,a] = butter(order, fs, band);
                butterEMG = filtfilt(b, a, EMGoutput);
                filtEMG = butterEMG - mean(butterEMG(1:fsamp));
        end
    end
end

end

