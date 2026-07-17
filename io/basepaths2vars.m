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
addParameter(p, 'basepaths', {});
addParameter(p, 'vars', string([]));
addParameter(p, 'flgPrnt', false, @islogical);
parse(p, varargin{:});

basepaths = p.Results.basepaths;
vars = p.Results.vars;
flgPrnt = p.Results.flgPrnt;

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
        if flgPrnt, warning('%s does not exist, skipping...', filepath), end
        continue
    end
    cd(filepath)

    for ifile = 1:length(vars)
        % Match <basename>.<var>.mat exactly, then fall back to
        % <basename>.<var>*.mat (e.g. a '.cellinfo' infix). The leading dot is
        % what keeps a loose substring from grabbing AccuSleep_states.mat when
        % vars{ifile} = 'sleep_states'. On a tie the shortest name wins.
        filename = dir(['*.', vars{ifile}, '.mat']);
        if isempty(filename)
            filename = dir(['*.', vars{ifile}, '*.mat']);
        end

        if isempty(filename)
            if flgPrnt, warning('No %s file in %s, skipping...', vars{ifile}, filepath), end
            v(ipath).(vars{ifile}) = [];
            continue
        elseif length(filename) > 1
            [~, iShort] = min(cellfun(@length, {filename.name}));
            if flgPrnt, warning('Multiple %s files; using %s.', vars{ifile}, filename(iShort).name); end
            filename = filename(iShort).name;
        else
            filename = filename(1).name;
        end

        % Load the file and dynamically assign the variable
        temp = load(filename);

        % Automatically detect variable name
        varNames = fieldnames(temp);
        if isscalar(varNames)
            v(ipath).(varNames{1}) = temp.(varNames{1});
        elseif any(strcmp(varNames, vars{ifile}))
            v(ipath).(vars{ifile}) = temp.(vars{ifile});
        else
            % Default to the first variable if multiple exist
            if flgPrnt, warning('Multiple variables in %s. Using first one.', filename), end
            v(ipath).(vars{ifile}) = temp.(varNames{1});
        end
    end
end

end

% EOF