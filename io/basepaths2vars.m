function v = basepaths2vars(varargin)
% Loads specified variables from multiple directories and organizes them in a struct array.
% 
% INPUT:
%   basepaths - cell array of directory paths.
%   vars      - string array of .mat file names (without extensions) to load.
% 
% OUTPUT:
%   v - struct array with fields corresponding to variables inside the .mat files.
%
% EXAMPLE:
%   v = basepaths2vars('basepaths', basepaths, 'vars', ["session"; "cell_metrics"]);
%
% DEPENDENCIES:
%   Requires that each basepath contains the specified .mat files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepaths', {}, @(x) iscell(x) && all(cellfun(@ischar, x)));
addParameter(p, 'vars', string([]), @isstring);
parse(p, varargin{:});

basepaths = p.Results.basepaths;
vars = p.Results.vars;

npaths = length(basepaths);

if isempty(basepaths) || isempty(vars)
    error('basepaths and vars must be specified.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Files from Basepaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear v

for ipath = 1 : npaths
    filepath = basepaths{ipath};
    if ~exist(filepath, 'dir')
        warning('%s does not exist, skipping...', filepath)
        continue
    end
    cd(filepath)
    [~, basename] = fileparts(filepath);
    
    for ifile = 1:length(vars)
        filename = dir(['*', vars{ifile}, '*.mat']);
        
        if length(filename) > 1
            warning('Multiple files with the name %s. Using first one %s.',...
                vars{ifile}, filename(1).name);
            filename = filename(1).name;
        elseif isempty(filename)
            warning('No %s file in %s, skipping...', vars{ifile}, filepath)
            v(ipath).(vars{ifile}) = [];
            continue
        else
            filename = filename(1).name;
        end
        
        % Load the file and dynamically assign the variable
        temp = load(filename);
        
        % Automatically detect variable name
        varNames = fieldnames(temp);
        if length(varNames) == 1
            v(ipath).(varNames{1}) = temp.(varNames{1});
        elseif any(strcmp(varNames, vars{ifile}))
            v(ipath).(vars{ifile}) = temp.(vars{ifile});
        else
            % Default to the first variable if multiple exist
            warning('Multiple variables in %s. Using first one.', filename)
            v(ipath).(vars{ifile}) = temp.(varNames{1});
        end
    end
end

end

% EOF