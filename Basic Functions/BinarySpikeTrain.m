function [binarySpT, smoothSpT] = BinarySpikeTrain(MotorUnits, AnalogSignal, varargin)
%This function takes the discahrge times of motor units and creates a
%binary spike train and a smoothed binary spike train

%   If trial is a ramo contraction, the smoothed binary spike train needs
%   to be the same length as an analog signal

if nargin > 1
    for i=1:2:size(varargin,2)
        switch varargin{i}
            case 'AnalogSignal'
                AnalogSignal = varargin{i+1};
            case 'highpass'
                HpHz = varargin{i+1};
            case 'windowLength'
                winL = varargin{i+1}./1000;
            otherwise
                errordlg('unknown argument to allChanEMG')
        end
    end
else
end

% AnalogSignal = [];
fsamp = 2048;
smoothSpT = [];

if ~exist('winL', 'var')
    winL = 1;
    mxWin = 2;
else 
    mxWin = winL./100;
end 

if ~exist('AnalogSignal', 'var') 
    for MU = 1:size(MotorUnits,2)
        index = MotorUnits{MU};
        binarySpT(MU,index) = 1;
    end

else
    
    trialLength = length(AnalogSignal);
    for MU = 1:size(MotorUnits,2)
        index = MotorUnits{MU};
        binarySpT(MU,index) = 1;
    end
    
    binarySpT(:,length(binarySpT):trialLength)=0; %adds a second of 0s to the end
    try
        smoothSpT= mxWin*fftfilt(hanning(round(fsamp*winL)),binarySpT')';
    catch
        disp('Could not create binary spike train');
    end
end

if exist('HpHz', 'var')
    smoothSpT = highpass(smoothSpT, HpHz,fsamp);
end 



end

