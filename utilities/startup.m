% startup.m
%
% script that launches automatically upon startup and whenever called upon
%
% place this file in userpath, e.g:
% D:\Program Files\Matlab\R2018b\toolbox\local\userpath.m
%
% 18 nov 19 LH

if ~isdeployed
    matlabpath(pathdef)
end

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

% EOF
