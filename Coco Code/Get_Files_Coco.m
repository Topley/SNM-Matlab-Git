function [filenames] = Get_Files_Coco(masterSheet, varargin)
%getCocoFiles Find trials given certain parameteres
%   subject = identify specific subject if only looking at single subjects
% removeSub = 1 = look at all subjects EXCEPT the identified subject, 0 =
% look at only the single identified subject
% rampNum = 0 get at all ramps, otherwise specify ramp of interest
% feedback must be specified as Ramp, Hold, Iramp
% muscle need to specify which muscle to files for
% Contraction - DF, PF, Coco

[num,txt,raw] = xlsread(masterSheet) ;
filenames = cellfun(@num2str,raw,'un',0);

subject = 'Coco99';
removeSub = 1;
rampNum = 0;
if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'subject'
        subject = varargin{i+1};
      case 'muscle'
        muscle = varargin{i+1};
      case 'feedback'
        fb = varargin{i+1};
      case 'rampNum'
        rampNum = varargin{i+1};
      case 'remove'
        removeSub = varargin{i+1};
      case 'contraction'
        contraction = varargin{i+1};
      otherwise
        errordlg('unknown argument')
    end
  end
else
end

if removeSub == 1
Subject2Remove = contains(filenames(:,1), subject);
filenames(Subject2Remove,:) = [];
else
Subject2Remove = ~contains(filenames(:,1), subject);
filenames(Subject2Remove,:) = [];
end 

muscle2Remove = ~contains(filenames(:,1), muscle) ;
filenames(muscle2Remove,:) = [];

feedback2Remove = ~contains(filenames(:,2), fb, 'IgnoreCase', false);
filenames(feedback2Remove,:) = [];

contraction2Remove = ~contains(filenames(:,3), contraction);
filenames(contraction2Remove,:) = [];

removeBadTrials = contains(filenames(:,5), 'Bad');
filenames(removeBadTrials,:) = [];

filenames = filenames(:,1);

end

