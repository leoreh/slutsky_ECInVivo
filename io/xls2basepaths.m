function basepaths = xls2basepaths(varargin)

% The XLS2BASEPATHS reads a session list from an Excel sheet, filters sessions
% based on positive and negative conditions (columns in the sheet), and
% constructs the full paths.
%
% INPUT (Optional Key-Value Pairs):
%   mname       (char) Mouse name which specifies which tab in the xls sheet
%               to load.
%   xlsname     (char) Full path to xls file with list of sessions.
%               Defaults to hardcoded path if empty.
%   pcond       (string array) Column names of logical values for each dir.
%               Only if true (1) will the session be loaded.
%               {"tempflag"}
%   ncond       (string array) Same as pcond but imposes a negative condition.
%               {"fepsp"}
%   dirCol      (string) Column name in xls sheet where dirnames exist.
%               "Session"
%   pathCol     (string) Column name in xls sheet where the path for each
%               dir exists.
%               "Mousepath"
%
% OUTPUT:
%   basepaths   (cell array) List of full directory paths.
%
% HISTORY:
%   21 oct 20 LH    updates
%   27 dec 21 LH    cell of structs instead of cell array (legacy)
%   10 dec 25 Agent Refactored from getSessionVars to xls2basepaths

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'mname', '', @ischar);
addOptional(p, 'xlsname', '', @ischar);
addOptional(p, 'pcond', "tempflag", @(x) isstring(x) || ischar(x) || iscell(x));
addOptional(p, 'ncond', "fepsp", @(x) isstring(x) || ischar(x) || iscell(x));
addOptional(p, 'dirCol', "Session", @(x) isstring(x) || ischar(x));
addOptional(p, 'pathCol', "Mousepath", @(x) isstring(x) || ischar(x));

parse(p, varargin{:})
mname = p.Results.mname;
xlsname = p.Results.xlsname;
pcond = string(p.Results.pcond);
ncond = string(p.Results.ncond);
dirCol = string(p.Results.dirCol);
pathCol = string(p.Results.pathCol);

%% ========================================================================
%  DEFAULTS
%  ========================================================================

% Default path to xls file with session metadata
if isempty(xlsname)
    xlsname = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Data summaries\sessionList.xlsx';
end

if ~isempty(xlsname)
    if ~ischar(xlsname) && ~contains(xlsname, 'xlsx')
        error('check xlsname')
    end
end

if isempty(mname)
    warning('Mouse name (mname) not specified. This might cause readtable to fail if sheet is not specified.');
end

%% ========================================================================
%  LOAD AND FILTER
%  ========================================================================

% Load xls
xlsfile = dir(xlsname);
if isempty(xlsfile)
    error('Excel file not found at: %s', xlsname);
end
sessionInfo = readtable(fullfile(xlsfile.folder, xlsfile.name), 'sheet', mname);
nrows = size(sessionInfo, 1);

% Select rows based on conditions
irow = true(nrows, 1);

% Positive conditions
for idir = 1 : length(pcond)
    icol = strcmp(sessionInfo.Properties.VariableNames, pcond(idir));
    if any(icol)
        irow = irow & sessionInfo{:, icol} == 1;
    end
end

% Negative conditions
for idir = 1 : length(ncond)
    icol = strcmp(sessionInfo.Properties.VariableNames, ncond(idir));
    if any(icol)
        irow = irow & sessionInfo{:, icol} ~= 1;
    end
end

%% ========================================================================
%  CONSTRUCT PATHS
%  ========================================================================

% Organize basepaths from dirnames and dirpaths
idirCol = strcmp(sessionInfo.Properties.VariableNames, dirCol);
ipathCol = strcmp(sessionInfo.Properties.VariableNames, pathCol);

dirnames = string(table2cell(sessionInfo(irow, idirCol)));
dirpaths = string(table2cell(sessionInfo(irow, ipathCol)));

% Check if we found any sessions
if isempty(dirnames)
    warning('No sessions found matching conditions for mouse %s', mname);
    basepaths = {};
    return;
end

% Construct full paths
basepaths = cell(1, length(dirnames));
% Assuming the first path can be used as default if individual paths are missing
default_path = dirpaths{1};

for ipath = 1 : length(dirnames)
    if isempty(dirpaths{ipath}) || ismissing(dirpaths{ipath})
        basepaths{ipath} = fullfile(default_path, dirnames{ipath});
    else
        basepaths{ipath} = fullfile(dirpaths{ipath}, dirnames{ipath});
    end

    % Convert to char for compatibility
    if isa(basepaths{ipath}, 'string')
        basepaths{ipath} = char(basepaths{ipath});
    end
end

end

% EOF
