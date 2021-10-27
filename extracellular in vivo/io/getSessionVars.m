function [varArray, dirnames, mousepath] = getSessionVars(varargin)

% loads specific matlab variables from specified directories and organizes
% them in a cell array. the directories can be specified directly as input
% or recieved from an xls file with certain conditions. assumes all
% variables in each dir are named "dirname.varname.mat".

% INPUT:
%   mousepath   string. path to mouse folder with all sessions {pwd}. if
%               empty will be retreived from xls
%   xlsname     string. name of xls file with list of sessions. must
%               include extension {'sessionList.xlsx'}
%   dirColumn 	string. column name in xls sheet where dirnames exist
%   vars        string array of matlab variables to load.
%   pcond       string. column name of logical values for each
%               session. only if true than session will be loaded. can be a
%               string array and than all conditions must be met.
%   ncond       string. same as pcond but imposes a negative condition.
%   sortDir     logical. sort varArray according to dirnames {true}
%
% DEPENDENCIES
%   none
%
% TO DO LIST:
%   add option to input dirnames
%
% 21 oct 20 LH      updates:
%   added default values for my convinience

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'mousepath', pwd);
addOptional(p, 'xlsname', '', @ischar);
addOptional(p, 'mname', '', @ischar);
addOptional(p, 'dirColumn', "Session", @isstring);
addOptional(p, 'vars', string([]), @isstring);
addOptional(p, 'pcond', "tempflag", @isstring);
addOptional(p, 'ncond', "fepsp", @isstring);
addOptional(p, 'sortDir', true, @islogical);
addOptional(p, 'dirnames', []);

parse(p, varargin{:})
mousepath       = p.Results.mousepath;
xlsname         = p.Results.xlsname;
mname           = p.Results.mname;
dirColumn       = p.Results.dirColumn;
vars            = p.Results.vars;
pcond           = p.Results.pcond;
ncond           = p.Results.ncond;
sortDir         = p.Results.sortDir;
dirnames        = p.Results.dirnames;

% default path to xls file with session metadata
if isempty(xlsname)
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
end

% deafuly variables to load. note this is congruent with assignVars.m
if isempty(vars)
    vars = ["session.mat";...
        "cell_metrics.cellinfo";...
        "spikes.cellinfo";...
        "fr.mat";...
        "datInfo";...
        "AccuSleep_states";...
        "sr.mat"];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get directory paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('dirnames', 'var') && isstring(dirnames)
    % ALT 1: user input dirnames
    dirnames = dirnames;
    
    % ALT 2: get dirnames from xlsx file
elseif ischar(xlsname) && contains(xlsname, 'xlsx')
    xlsfile = dir(xlsname);
    sessionInfo = readtable(fullfile(xlsfile.folder, xlsfile.name),...
        'sheet', mname);
    icol = strcmp(sessionInfo.Properties.VariableNames, dirColumn);
    dirnames = string(table2cell(sessionInfo(:, icol)));
    basepaths = string(table2cell(sessionInfo(:, 1)));
    basepath_default = basepaths{1};
    
    % check dirnames meet conditions
    clear irow icol
    irow = ones(length(dirnames), 1);
    for i = 1 : length(pcond)
        icol = strcmp(sessionInfo.Properties.VariableNames, char(pcond(i)));
        if any(icol)
            irow = irow & sessionInfo{:, icol} == 1;
        end
    end
    for i = 1 : length(ncond)
        icol = strcmp(sessionInfo.Properties.VariableNames, char(ncond(i)));
        if any(icol)
            irow = irow & sessionInfo{:, icol} ~= 1;
        end
    end
    dirnames = dirnames(irow);
    dirnames(strlength(dirnames) < 1) = [];
    basepaths = basepaths(irow);
end

if sortDir
    dirnames = string(natsort(cellstr(dirnames)));
end
ndirs = length(dirnames);

% load files
varArray = cell(length(dirnames), length(vars));
for i = 1 : ndirs
    
    % handle cases where not all sessions are in the same basepath
    if exist('basepaths', 'var')
        if ~isempty(basepaths{i})
            mousepath = basepaths{i};
        else
            mousepath = basepath_default;
        end
    else
        mousepath = fileparts(pwd);
    end
    
    filepath = char(fullfile(mousepath, dirnames(i)));
    if ~exist(filepath, 'dir')
        warning('%s does not exist, skipping...', filepath)
        continue
    end
    cd(filepath)
    
    for ii = 1 : length(vars)
        filename = dir(['*', vars{ii}, '*']);
        if length(filename) == 1
            varArray{i, ii} = load(filename.name);
        elseif length(filename) > 1
            varArray{i, ii} = load(filename(2).name);       % lazy fix for tdt if both Raw and EMG streams were loaded
        else
            warning('no %s file in %s, skipping', vars{ii}, filepath)
        end
    end
end


end