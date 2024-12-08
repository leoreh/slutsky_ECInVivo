function [varArray, basepaths] = getSessionVars(varargin)

% loads a list of specific variables from multiple directories and
% organizes them in struct. the directories (basepaths) can be
% specified directly as input or recieved from an xls file with certain
% conditions. assumes all variables in each dir are named
% "dirname.varname.mat".

% INPUT:
%   basepaths   cell of chars to dirs of interest.
%   xlsname     char. name of xls file with list of sessions. must
%               include extension {'sessionList.xlsx'}
%   mname       char. mouse name which specifies which tab in the xls sheet
%               to load
%   dirCol      char. column name in xls sheet where dirnames exists
%   pathCol     char. column name in xls sheet where the path for each die
%               exists
%   varsFile    string array of matlab variables to load.
%   varsName    string array of the struct field names to place each var.
%               see example below
%   pcond       char. column name of logical values for each
%               dir. only if true than session will be loaded. can be a
%               string array and than all conditions must be met.
%   ncond       string. same as pcond but imposes a negative condition.
%
% EXAMPLE
%   basepaths = [{'G:\Data\lh93\lh93_210811_102035'},
%   {'G:\Data\lh93\lh93_210813_110609'};
%   varArray = getSessionVars('basepaths', basepaths,...
%       'varsFile', ["cell_metrics"; "spktimesMetrics"],...
%       'varsName', ["cm"; "st"])
%   The output varArray is a struct of length 2 (one for each basepath)
%   with fields .cm and .st that contain the variables loaded from
%   cell_metrics and spktimesMetrics, respectively.
% 
% DEPENDENCIES
%   none
%
% TO DO LIST:
%   add option to input dirnames (done)
%
% 21 oct 20 LH      updates:
% 27 dec 21 LH      cell of structs instead of cell array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepaths', []);
addOptional(p, 'xlsname', '', @ischar);
addOptional(p, 'mname', '', @ischar);
addOptional(p, 'dirCol', "Session", @isstring);
addOptional(p, 'pathCol', "Mousepath", @isstring);
addOptional(p, 'varsFile', string([]), @isstring);
addOptional(p, 'varsName', string([]), @isstring);
addOptional(p, 'pcond', "tempflag", @isstring);
addOptional(p, 'ncond', "fepsp", @isstring);

parse(p, varargin{:})
basepaths       = p.Results.basepaths;
xlsname         = p.Results.xlsname;
mname           = p.Results.mname;
dirCol          = p.Results.dirCol;
pathCol         = p.Results.pathCol;
varsFile        = p.Results.varsFile;
varsName        = p.Results.varsName;
pcond           = p.Results.pcond;
ncond           = p.Results.ncond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default path to xls file with session metadata
if isempty(xlsname)
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
end

% deafult variables to load
if isempty(varsFile)
    varsFile = ["session"; "cell_metrics.cellinfo"; "spikes.cellinfo";...
        "fr"; "datInfo"; "AccuSleep_states"; "sr"; "st_metrics";
        "swv_metrics"; "ripp"];
end

if isempty(varsName)
    varsName = ["session"; "cm"; "spikes"; "fr"; "datInfo"; "ss";...
        "sr"; "st"; "swv"; "ripp"];
end

if length(varsFile) ~= length(varsName)
    error('varsFile and varsName must be of equal length')
end

if ~isempty(xlsname)
    if ~ischar(xlsname) && ~contains(xlsname, 'xlsx')
        error('check xlsname')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get directory paths from xls 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(basepaths)
    
    % load xls
    xlsfile = dir(xlsname);
    sessionInfo = readtable(fullfile(xlsfile.folder, xlsfile.name),...
        'sheet', mname);
    nrows = size(sessionInfo, 1);

    % select only rows in the xls that meet the positive and negative
    % conditions conditions
    clear irow icol
    irow = ones(nrows, 1);
    for idir = 1 : length(pcond)
        icol = strcmp(sessionInfo.Properties.VariableNames, char(pcond(idir)));
        if any(icol)
            irow = irow & sessionInfo{:, icol} == 1;
        end
    end
    for idir = 1 : length(ncond)
        icol = strcmp(sessionInfo.Properties.VariableNames, char(ncond(idir)));
        if any(icol)
            irow = irow & sessionInfo{:, icol} ~= 1;
        end
    end   
    
    % organize basepaths from dirnames and dirpaths
    idirCol = strcmp(sessionInfo.Properties.VariableNames, dirCol);
    ipathCol = strcmp(sessionInfo.Properties.VariableNames, pathCol);  
    dirnames = string(table2cell(sessionInfo(irow, idirCol)));
    dirpaths = string(table2cell(sessionInfo(irow, ipathCol)));
    default_path = dirpaths{1};   
    for ipath = 1 : length(dirnames)
        if isempty(dirpaths{ipath})
            basepaths{ipath} = fullfile(default_path, dirnames{ipath});
        else
            basepaths{ipath} = fullfile(dirpaths{ipath}, dirnames{ipath});
        end
    end
end
ndirs = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files and organize in cell struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load files
clear varArray
for idir = 1 : ndirs
    
    filepath = basepaths{idir};
    if ~exist(filepath, 'dir')
        warning('%s does not exist, skipping...', filepath)
        continue
    end
    cd(filepath)
    [~, basename] = fileparts(filepath);
    
    for ifile = 1 : length(varsFile)
        filename = dir(['*', varsFile{ifile}, '*']);
        if length(filename) > 1
            % if more than one file share the same name, default is
            % to load the first. specific corrections can be applied below
            if contains(varsFile{ifile}, 'datInfo')             % tdt has a datInfo for each stream
                fileidx = contains({filename.name}, 'EMG');
                filename = filename(~fileidx).name;
            elseif any(~contains({filename.name}, basename))     % mea raw data names is a mess
                filename = filename(contains({filename.name}, basename)).name;
            elseif any(~contains({filename.name}, ['.', varsFile{ifile}, '.mat']))     % e.g. fr instead of "from" 
                filename = filename(contains({filename.name}, ['.', varsFile{ifile}, '.mat'])).name;
            else
                filename = filename(1).name;
            end
        elseif isempty(filename)         
            warning('no %s file in %s, skipping', varsFile{ifile}, filepath)
            varArray(idir).(varsName{ifile}) = [];
            continue
        else
            filename = filename(1).name;
        end
        
        % correct special cases where var name does not fit file name or
        % field name
        temp = load(filename);
        if strcmp(varsName{ifile}, 'cm')
                varArray(idir).(varsName{ifile}) = temp.cell_metrics;
        elseif strcmp(varsName{ifile}, 'sr')
            varArray(idir).(varsName{ifile}) = temp.fr;
        elseif strcmp(varsName{ifile}, 'psdEmg')
            varArray(idir).(varsName{ifile}) = temp.psd;
        else
            varArray(idir).(varsName{ifile}) = temp.(varsName{ifile});
        end
    end
end

end