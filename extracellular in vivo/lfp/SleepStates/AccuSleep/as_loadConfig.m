function [cfg_colors, cfg_names, cfg_weights] = as_loadConfig(configfile)

% loads accuSleep params from configuration file. full name nad path of
% config file can be specified, if not will look in default location, if
% not will search in same path as script, if not will ask user to find file
% example call: [~, cfg_names, ~] = as_loadConfig([]);

% 08 jun 21 LH

if isempty(configfile)
    configfile = 'D:\Code\slutskycode\extracellular in vivo\lfp\SleepStates\AccuSleep\AS_config.mat';
end
if ~exist(configfile, 'file')
    scriptfile = mfilename('fullpath');
    scriptpath = fileparts(scriptfile);
    configfile = fullfile(scriptpath, 'AS_config.mat');
    
    if ~exist(configfile, 'file')
        [configfile, configpath] = uigetfile('', 'Please select  the configuration file');
        configfile = [fullfile(configpath, configfile), '.mat'];
    end
end
load(configfile)

% assignin('base', 'cfg_colors', cfg_colors)
% assignin('base', 'cfg_names', cfg_names)
% assignin('base', 'cfg_weights', cfg_weights)

end