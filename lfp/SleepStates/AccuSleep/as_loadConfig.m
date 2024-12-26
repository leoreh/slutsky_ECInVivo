function [cfg] = as_loadConfig(varargin)

% loads accuSleep params from configuration file. full name and path of
% config file can be specified. if not will look in default location, if
% not will search in same path as script, if not will ask user to find file
% example call: cfg = as_loadConfig();

% 08 jun 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'configfile', []);
addOptional(p, 'flgEmg', false, @islogical);

parse(p, varargin{:})
configfile      = p.Results.configfile;
flgEmg          = p.Results.flgEmg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(configfile)
    if flgEmg
        configfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\as_configEmg.mat';
    else
        configfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\as_config.mat';
    end
end
if ~exist(configfile, 'file')
    scriptfile = mfilename('fullpath');
    scriptpath = fileparts(scriptfile);
    configfile = fullfile(scriptpath, 'as_config.mat');
    
    if ~exist(configfile, 'file')
        [configfile, configpath] = uigetfile('', 'Please select  the configuration file');
        configfile = [fullfile(configpath, configfile)];
    end
end
load(configfile)

end