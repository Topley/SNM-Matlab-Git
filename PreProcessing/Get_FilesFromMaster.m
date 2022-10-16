function [filenames] = Get_FilesFromMaster(masterSheet, varargin)
%%% Get_FilesFromMaster Find trials given the parameteres below
%   subject = Specify the subject if to get only that subjects files 
%   removeSub = Input 1 to remove a specific subject from the files of interest 
%   trialNumber = Get a specific trial when there are multiples of a trial
%   trial = Specify the specific trial type if interested in a specific kind of trial

[num,txt,raw] = xlsread(masterSheet);
filenames = cellfun(@num2str,raw,'un',0);

if nargin
    for i=1:2:size(varargin,2)
        switch varargin{i}
            case 'subject'
                subject = varargin{i+1};
            case 'trialType'
                trialType = varargin{i+1};
            case 'trialNumber'
                trialNumber = varargin{i+1};
            case 'remove'
                removeSub = varargin{i+1};
            otherwise
                errordlg('unknown argument')
        end
    end
else
end

if exist('removeSub', 'var') == 1
    Subject2Remove = contains(filenames(:,2), subject);
    filenames(Subject2Remove,:) = [];
else
    try
        Subject2Remove = ~contains(filenames(:,2), subject);
        filenames(Subject2Remove,:) = [];
    end
end

if exist('trialType', 'var') == 1
    type2Remove = ~contains(filenames(:,4), trialType, 'IgnoreCase', false);
    filenames(type2Remove,:) = [];
end

if exist('trialNumber', 'var') == 1
    number2Remove = ~contains(filenames(:,5), trialNumber, 'IgnoreCase', false);
    filenames(number2Remove,:) = [];
end

removeBadTrials = contains(filenames(:,6), 'Bad'); % remove files that have "bad" in the comments indicating an issue during collection
filenames(removeBadTrials,:) = [];

filenames = filenames(:,[2:4]); % returns the associated c3d file in column 1, the OTB file in column 2, and the trial type in column 3

end
