% startup.m
% 
% script that launches automatically upon startup and whenever called upon
% 
% place this file in userpath, e.g:
% D:\Program Files\Matlab\R2018b\toolbox\local\userpath.m)
% 
% 18 nov 19 LH

matlabpath(pathdef)

% command line formatting
format compact
format long g
dbstop if error

% clear some variables
diary('off')
close all
clear all
clc

% EOF
