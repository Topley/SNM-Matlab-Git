function [] = startup()
%Function to add paths so the startup file does not need to be added to
%scripts

%   add tableau graphics to filepath
addpath(genpath('C:\SNM Matlab SDK'));

% Add folders containing matlab code
% addpath(genpath(uigetdir('C:\Users\tuc43377\Desktop\CoCo project\Matlab Code')));
tableau = Tableau(); %loading tableau colormap


% set default settings for figures etc
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultFigureColorMap', tableau);
set(0, 'DefaultAxesColorOrder', tableau);
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');

pos = [1775 485 560 420];
set(0, 'DefaultFigurePosition', pos);

% folder to begin in at startup
%cd(uigetdir('C:\Users\tuc43377\Desktop\CoCo project'))
cd('G:\TRD Hub');


end

