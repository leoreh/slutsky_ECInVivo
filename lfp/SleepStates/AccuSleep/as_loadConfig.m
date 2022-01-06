function [cfg] = as_loadConfig(configfile)

% loads accuSleep params from configuration file. full name and path of
% config file can be specified. if not will look in default location, if
% not will search in same path as script, if not will ask user to find file
% example call: cfg = as_loadConfig();

% 08 jun 21 LH

if nargin == 0
    configfile = [];
end

if isempty(configfile)
    configfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\as_config.mat';
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