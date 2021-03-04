% startup.m
%
% script that launches automatically upon startup and whenever called upon
%
% place this file in userpath, e.g:
% D:\Program Files\Matlab\R2018b\toolbox\local\userpath.m
%
% 18 nov 19 LH

matlabpath(pathdef)

% command line formatting
format compact
format long g
dbstop if error

% clear workspace and command window
diary('off')
close all
clear all
close('all')
clc

% exist debugging mode
if feature('IsDebugMode')
    dbquit all
end

% Get rid of warning about directory already existing:
% "Warning: Directory already exists."
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% change default graphics
% NOTE groot is the root of all graphical objects (previously referred to
% as 0)
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');
set(groot, 'DefaultLegendInterpreter', 'none');
% set(groot, 'DefaultAxesFontName', 'Diverda Sans Com')
% set(groot, 'DefaultTextFontName', 'Diverda Sans Com')
% set(groot, 'DefaultAxesFontSize', 12)
% set(groot, 'DefaultAxesOuterPosition', [0 0 1 1])
% set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.2)
% set(groot, 'DefaultAxesFontWeight', 'normal')
% set(groot, 'DefaultAxesTitleFontWeight', 'bold')
% set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.4)
% 
% 
% set(groot, 'DefaultFigurePaperUnits', 'normalize')
% set(groot, 'DefaultFigureUnits', 'normalize')
% set(groot, 'DefaultFigurePaperPosition', [0 0 1 1])
% set(groot, 'DefaultFigurePosition', [0.2 0.075 0.6 0.85])
% set(groot, 'DefaultFigurePaperType', 'a4')
% set(groot, 'DefaultFigureColor', [1 1 1])       % figure white instead of gray
% set(groot, 'DefaultFigurePaperOrientation', 'portrait')
% 
% set(groot, 'DefaultLineLineWidth', 0.01)


% factory values
% set(groot, 'DefaultTextInterpreter', 'tex');
% set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
% set(groot, 'DefaultLegendInterpreter', 'tex');
% set(groot, 'DefaultAxesFontName', 'Helvetica')
% set(groot, 'DefaultTextFontName', 'Helvetica')
% set(groot, 'DefaultAxesFontSize', 10)
% set(groot, 'DefaultAxesOuterPosition', [0 0 1 1])
% set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1)
% set(groot, 'DefaultAxesFontWeight', 'normal')
% set(groot, 'DefaultAxesTitleFontWeight', 'bold')
% set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.1)
% 
% 
% set(groot, 'DefaultFigurePaperUnits', 'inches')
% set(groot, 'DefaultFigureUnits', 'pixels')
% set(groot, 'DefaultFigurePaperPosition', [0.25 2.5 8 6])
% set(groot, 'DefaultFigurePosition', [200 200 600 500])
% set(groot, 'DefaultFigurePaperType', 'usletter')
% set(groot, 'DefaultFigureColor', [0.94 0.94 0.94])       % figure white instead of gray
% set(groot, 'DefaultFigurePaperOrientation', 'portrait')
% 
% set(groot, 'DefaultLineLineWidth', 0.5)






% EOF
