function [] = startup()
%Function to add paths so the startup file does not need to be added to
%scripts

%   add tableau graphics to filepath
addpath(genpath('C:\SNM Git Project'));

% Add folders containing matlab code
% addpath(genpath(uigetdir('C:\Users\tuc43377\Desktop\CoCo project\Matlab Code')));
tableau = Tableau(); %loading tableau colormap


% set default settings for figures etc
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultFigureColorMap', tableau);
set(0, 'DefaultAxesColorOrder', tableau);
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');

set(0,'units','pixels');
screensize = get(0,'ScreenSize');

pos = [screensize(3).*.75,   30, screensize(3).*.25,  screensize(3).*.25];
set(0, 'DefaultFigurePosition', pos);

% folder to begin in at startup
cd(uigetdir('C:'))
rmpath('C:\SNM Git Project\zOld Matlab Code')


end

