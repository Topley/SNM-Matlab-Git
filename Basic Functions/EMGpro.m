function [varargout] = EMGpro(arrayEMG, varargin)
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
 
if size(varargin,2) > 0
    for i=1:2:size(varargin,2)
        switch varargin{i} % parsing inputs
            case 'channel'
                channel = varargin{i+1};
            case 'filter'
                filtVars = varargin{i+1};
                filtType = filtVars{1};
                switch filtType
                    case 'rawRMS'
                        rmsL = filtVars{2}; %unfiltered RMS length 
                    case 'filtRMS' 
                        rmsL = filtVars{2}; %RMS length
                        rmsWindow = filtVars{3}; %filter window length
                    case 'butter'
                        order = filtVars{1}; % filter order
                        fs = filtVars{2}; % frequency cutoff/band
                        band = filtVars{3}; % filter type - high, low, bandpass
                end
                
            otherwise
                errordlg('unknown argument to allChanEMG')
        end
    end
end

if exist('channel', 'var')
    EMGoutput = abs(arrayEMG(channel, :));
    varargout{1} = EMGoutput;
else
    EMGoutput = (abs(sum(arrayEMG))) ./size(arrayEMG, 1);
    varargout{1} = EMGoutput;
end

try
    switch filtType
        case 'rawRMS'
            varargout{2} = sqrt(rms(EMGoutput .^2, rmsL));
        case 'filtRMS'
            varargout{2} = sqrt(movmean(rms(EMGoutput .^2, rmsL), rmsWindow));
        case 'butter'
            [b,a] = butter(order, fs, band);
            varargout{2} = filtfilt(b, a, varargout{1});
    end
catch
    disp('Error in EMGpro Function')
end

end

