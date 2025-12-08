function basepaths = xls2basepaths(varargin)

% Extracts directory paths (basepaths) from an Excel sheet based on conditions.
% 
% INPUT:
%   xlsname  - char. Excel file name with session list, including extension.
%   mname    - char. Mouse name (sheet to load from the Excel file).
%   dirCol   - char. Column name for directory names.
%   pathCol  - char. Column name for paths.
%   pcond    - string array. Positive condition column names (must be true).
%   ncond    - string array. Negative condition column names (must be false).
% 
% OUTPUT:
%   basepaths - cell array of full paths.
%
% EXAMPLE:
%   basepaths = xls2Basepaths('mname', 'lh132', 'pcond', ["tempflag"]);
%
% DEPENDENCIES:
%   Requires Excel file with appropriate format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'xlsname', 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Data summaries\sessionList.xlsx', @ischar);
addParameter(p, 'mname', '', @ischar);
addParameter(p, 'dirCol', "Session", @isstring);
addParameter(p, 'pathCol', "Mousepath", @isstring);
addParameter(p, 'pcond', "", @isstring);
addParameter(p, 'ncond', "", @isstring);

parse(p, varargin{:});
xlsname     = p.Results.xlsname;
mname       = p.Results.mname;
dirCol      = p.Results.dirCol;
pathCol     = p.Results.pathCol;
pcond       = p.Results.pcond;
ncond       = p.Results.ncond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Excel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlsfile = dir(xlsname);
sessionInfo = readtable(fullfile(xlsfile.folder, xlsfile.name), 'Sheet', mname);
nrows = size(sessionInfo, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irow = ones(nrows, 1);

for idir = 1:length(pcond)
    icol = strcmp(sessionInfo.Properties.VariableNames, char(pcond(idir)));
    if any(icol)
        irow = irow & sessionInfo{:, icol} == 1;
    end
end

for idir = 1:length(ncond)
    icol = strcmp(sessionInfo.Properties.VariableNames, char(ncond(idir)));
    if any(icol)
        irow = irow & sessionInfo{:, icol} ~= 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Basepaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirnames = string(table2cell(sessionInfo(irow, strcmp(sessionInfo.Properties.VariableNames, dirCol))));
dirpaths = string(table2cell(sessionInfo(irow, strcmp(sessionInfo.Properties.VariableNames, pathCol))));

default_path = dirpaths{1};

for ipath = 1:length(dirnames)
    if isempty(dirpaths{ipath})
        basepaths{ipath} = fullfile(default_path, dirnames{ipath});
    else
        basepaths{ipath} = fullfile(dirpaths{ipath}, dirnames{ipath});
    end
end

end
