function [v, tGlobal] = mea_tAlign(v, varargin)
% MEA_TALIGN Aligns firing rate matrices and other variables to perturbation onset.
%
%   [v, t] = MEA_TALIGN(v, ...) receives a struct array 'v' (output of
%   basepaths2vars) and aligns specified variables such that the
%   perturbation index (idxPert) is identical across all files. It also
%   returns a global time vector 't'.
%
%   INPUTS:
%       v           - (struct array) Data structure containing variables to align.
%
%   OPTIONAL (Key-Value Pairs):
%       varMap      - (struct) Mapping of variable names to align.
%                     Fields are arbitrary names, values are paths in 'v'.
%                     Example: varMap.fr = 'mea.fr'; varMap.dyn = 'dyn';
%                     Default: assumes 'fr' exists in v.
%       refVar      - (char) Key in varMap to use as reference for idxPert.
%                     Default: 'fr'.
%
%   OUTPUTS:
%       v           - (struct array) Updated struct with aligned fields.
%       tGlobal     - (vector) Global time vector [seconds].
%
%   See also: MVP_FRPREP, BASEPATHS2VARS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'v', @isstruct);
addParameter(p, 'varMap', struct('fr', 'fr'), @isstruct);
addParameter(p, 'refVar', 'fr', @ischar);

parse(p, v, varargin{:});
varMap = p.Results.varMap;
refVar = p.Results.refVar;


%% ========================================================================
%  ANALYZE ALIGNMENT
%  ========================================================================

maxPre = 0;
maxPost = 0;
binSize = [];
refPath = varMap.(refVar);

% Determine the maximum pre- and post-perturbation durations
for iFile = 1:length(v)

    % Access reference variable
    % We need to find the 'info' struct. If refPath points to a field inside
    % a struct (e.g. 'fr.fr'), we assume 'info' is a sibling (e.g. 'fr.info').
    % If refPath is the struct itself (e.g. 'fr'), then 'info' is a child.

    try
        [val, parent] = getFieldAndParent(v(iFile), refPath);
    catch
        continue;
    end

    if isempty(val), continue; end
    iInfo = [];
    nBins = 0;

    % Strategies to find info
    if isstruct(val) && isfield(val, 'info')
        iInfo = val.info;
        if isfield(val, 'fr')
            nBins = size(val.fr, 2);
        else
            % Fallback if no .fr field but is struct
            warning('Reference variable is a struct but lacks .fr field.');
            % Try to guess size from largest numeric field? no, unsafe.
            continue;
        end
    elseif isstruct(parent) && isfield(parent, 'info')
        iInfo = parent.info;
        % If var is specific (e.g. fr.fr), use its size
        nBins = size(val, 2);
    end

    if isempty(iInfo)
        warning('Could not find info struct for file %d via reference %s', iFile, refVar);
        continue;
    end

    % Check consistency of binSize
    if isempty(binSize)
        binSize = iInfo.binSize;
    elseif iInfo.binSize ~= binSize
        warning('Bin size mismatch in file %d. Expected %f, got %f.', iFile, binSize, iInfo.binSize);
    end

    idxPert = iInfo.idxPert;

    nPre = idxPert - 1;
    nPost = nBins - idxPert;

    if nPre > maxPre
        maxPre = nPre;
    end

    if nPost > maxPost
        maxPost = nPost;
    end
end

% Calculate global dimensions
globalNBins = maxPre + 1 + maxPost;
globalIdxPert = maxPre + 1;


%% ========================================================================
%  CREATE GLOBAL TIME VECTOR
%  ========================================================================

if isempty(binSize)
    tGlobal = [];
    warning('No valid reference data found.');
    return;
end

% Create 0-centered index vector
tIdx = (1:globalNBins) - globalIdxPert;

% Convert to initial physical time (seconds)
tSec = tIdx * binSize;

% Apply corrections (matching logic in mea_frPrep)
tGlobal = zeros(size(tSec));

% Pre-perturbation (x3 scaling)
tGlobal(tIdx < 0) = tSec(tIdx < 0) * 3;

% Post-perturbation (x6 scaling)
tGlobal(tIdx >= 0) = tSec(tIdx >= 0) * 6;


%% ========================================================================
%  ALIGN DATA
%  ========================================================================

varKeys = fieldnames(varMap);

for iFile = 1:length(v)

    % --- Get Alignment Info for this File ---
    try
        [val, parent] = getFieldAndParent(v(iFile), refPath);
    catch
        continue;
    end

    if isempty(val), continue; end

    % Find info again (same logic as above)
    iInfo = [];
    nBinsRef = 0;

    if isstruct(val) && isfield(val, 'info')
        iInfo = val.info;
        if isfield(val, 'fr')
            nBinsRef = size(val.fr, 2);
        end
    elseif isstruct(parent) && isfield(parent, 'info')
        iInfo = parent.info;
        nBinsRef = size(val, 2);
    end

    if isempty(iInfo) || nBinsRef == 0, continue; end

    idxPert = iInfo.idxPert;
    nPre = idxPert - 1;
    padPre = maxPre - nPre;
    currentPost = nBinsRef - idxPert;
    padPost = maxPost - currentPost;


    % --- Align All Variables ---
    for iKey = 1:length(varKeys)
        key = varKeys{iKey};
        path = varMap.(key);

        try
            data = getFieldFromPath(v(iFile), path);
        catch
            continue;
        end

        if isempty(data), continue; end

        % Logic based on data type
        if isnumeric(data)
            % If numeric, check dimensions against reference bins
            if size(data, 2) == nBinsRef
                nRows = size(data, 1);
                valAligned = [nan(nRows, padPre), data, nan(nRows, padPost)];
                v(iFile) = setFieldFromPath(v(iFile), path, valAligned);
            end

        elseif isstruct(data)
            % If struct, iterate fields (Legacy / Whole-Struct Mode)
            flds = fieldnames(data);
            for iFld = 1:length(flds)
                fldName = flds{iFld};
                val = data.(fldName);

                % Check name first to prioritize time vector replacement
                if strcmp(fldName, 't') || strcmp(fldName, 'time')
                    data.(fldName) = tGlobal;

                elseif isnumeric(val) && ~isempty(val) && size(val, 2) == nBinsRef
                    nRows = size(val, 1);
                    valAligned = [nan(nRows, padPre), val, nan(nRows, padPost)];
                    data.(fldName) = valAligned;

                elseif strcmp(fldName, 'info') && strcmp(key, refVar)
                    data.info.idxPert = globalIdxPert;
                    data.info.tAlign = true;
                end
            end

            % Update struct in v
            v(iFile) = setFieldFromPath(v(iFile), path, data);
        end
    end
end

end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function val = getFieldFromPath(s, path)
parts = strsplit(path, '.');
val = s;
for i = 1:length(parts)
    val = val.(parts{i});
end
end

function [val, parent] = getFieldAndParent(s, path)
% Returns value and its immediate parent
parts = strsplit(path, '.');
val = s;
parent = [];
for i = 1:length(parts)
    parent = val;
    val = val.(parts{i});
end
end

function s = setFieldFromPath(s, path, val)
parts = strsplit(path, '.');
s = setFieldRecursive(s, parts, val);
end

function s = setFieldRecursive(s, parts, val)
if length(parts) == 1
    s.(parts{1}) = val;
else
    s.(parts{1}) = setFieldRecursive(s.(parts{1}), parts(2:end), val);
end
end
