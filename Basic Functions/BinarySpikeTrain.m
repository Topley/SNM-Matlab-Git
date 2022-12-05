function [binSpikes, smoothBST] = BinarySpikeTrain(MotorUnits, varargin)
%This function takes the discahrge times of motor units and creates a
%binary spike train and a smoothed binary spike train

%   If trial is a ramo contraction, the smoothed binary spike train needs
%   to be the same length as an analog signal
fsamp = 2048;
smoothBST = [];
winL = 1;
mxWin = 2;

if nargin > 1
    for i=1:2:size(varargin,2)
        switch varargin{i}
            case 'reference'
                referenceSignal = varargin{i+1};
            case 'filter'
                filtParams = varargin{i+1};
                    winL = filtParams(1);
                    mxWin = filtParams(2);
            otherwise
                errordlg('unknown argument to allChanEMG')
        end
    end
else
end

if empty(referenceSignal)
    for MU = 1:size(MotorUnits,2)
        index = MotorUnits{MU};
        binSpikes(MU,index) = 1;
    end
    trialLength = size(binSpikes,2);
    smoothBST= mxWin*fftfilt(hanning(round(fsamp*winL)),binSpikes');
else
    
    trialLength = length(referenceSignal);
    for MU = 1:size(MotorUnits,2)
        index = MotorUnits{MU};
        binSpikes(MU,index) = 1;
    end
    binSpikes(:,length(binSpikes):trialLength)=0; %adds a second of 0s to the end
end        

if nargout > 1
smoothBST= mxWin*fftfilt(hanning(round(fsamp*winL)),binSpikes')';
binSpikes = binSpikes';
else
binSpikes = binSpikes';
end

end

