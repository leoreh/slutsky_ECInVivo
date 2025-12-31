function [v, tGlobal] = mea_tAlign(v, varMap, pertMap, varargin)
% MEA_TALIGN Aligns firing rate matrices and other variables to perturbation onset.
%
%   [v, t] = MEA_TALIGN(v, varMap, pertMap, ...) receives a struct array 'v'
%   and aligns specified variables such that the perturbation index (idxPert)
%   is identical across all files. It also returns a global time vector 't'.
%
%   INPUTS:
%       v           - (struct array) Data structure containing variables to align.
%       varMap      - (struct) Mapping of variable names to align.
%                     Fields are arbitrary names, values are paths in 'v'.
%                     Example: varMap.fr = 'mea.fr'; varMap.dyn = 'dyn';
%       pertMap     - (char) Path to the perturbation index variable in 'v'.
%                     Example: 'fr.info.idxPert'.
%       varargin    - (param/value) Optional parameters:
%                     'binSize' : (num) Time bin size in seconds {60}
%                     'flgEdge' : (log) Alignment mode {true}
%                                 true  : Intersection (clip to smallest window)
%                                 false : Union (pad to largest window)
%
%   OUTPUTS:
%       v           - (struct array) Updated struct with aligned fields.
%       tGlobal     - (vector) Global time vector [seconds].
%
%   See also: MEA_FRPREP, BASEPATHS2VARS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'v', @isstruct);
addRequired(p, 'varMap', @isstruct);
addRequired(p, 'pertMap', @ischar);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'flgEdge', true, @islogical);

parse(p, v, varMap, pertMap, varargin{:});
binSize = p.Results.binSize;
flgEdge = p.Results.flgEdge;

%% ========================================================================
%  ANALYZE ALIGNMENT DIMENSIONS
%  ========================================================================

maxPre = 0;
maxPost = 0;
minPre = inf;
minPost = inf;

% Determine the pre- and post-perturbation durations for all files
for iFile = 1:length(v)

    % Access perturbation index
    idxPert = getFieldFromPath(v(iFile), pertMap);

    if isempty(idxPert)
        warning('Could not find idxPert for file %d at %s', iFile, pertMap);
        continue;
    end

    % Determine nBins from the first valid numeric variable in varMap
    nBins = 0;
    varKeys = fieldnames(varMap);

    for iVar = 1:length(varKeys)
        path = varMap.(varKeys{iVar});
        val = getFieldFromPath(v(iFile), path);

        if isnumeric(val) && ~isempty(val)
            nBins = size(val, 2);
            break;
        elseif isstruct(val)
            % Try to find existing numeric fields inside struct (e.g. .fr)
            flds = fieldnames(val);
            for f = 1:length(flds)
                dd = val.(flds{f});
                if isnumeric(dd) && ~isempty(dd) && size(dd, 2) > 1
                    nBins = size(dd, 2);
                    break;
                end
            end
            if nBins > 0, break; end
        end
    end

    if nBins == 0
        warning('Could not determine nBins for file %d', iFile);
        continue;
    end

    % Logic:
    % Pre-pert bins: 1 to (idxPert - 1)
    % Post-pert bins: idxPert to nBins (including pert bin)

    currentPre = idxPert - 1;
    currentPost = nBins - idxPert + 1;

    maxPre = max(maxPre, currentPre);
    maxPost = max(maxPost, currentPost);
    minPre = min(minPre, currentPre);
    minPost = min(minPost, currentPost);
end

% Set Target Dimensions based on flgEdge
if flgEdge
    % INTERSECTION: Smallest common window
    targetPre = minPre;
    targetPost = minPost;
else
    % UNION: Largest window
    targetPre = maxPre;
    targetPost = maxPost;
end

% Global Total Bins
globalNBins = targetPre + targetPost;
% Global Index of Perturbation
globalIdxPert = targetPre + 1;


%% ========================================================================
%  CREATE GLOBAL TIME VECTOR
%  ========================================================================

% Create 0-centered index vector
tIdx = (1:globalNBins) - globalIdxPert;

% Convert to initial physical time (seconds)
tSec = tIdx * binSize;

% Apply corrections (matching logic in mea_frPrep)
tGlobal = zeros(size(tSec));

% Pre-perturbation (x3 scaling from legacy code)
tGlobal(tIdx < 0) = tSec(tIdx < 0) * 3;

% Post-perturbation (x6 scaling from legacy code)
tGlobal(tIdx >= 0) = tSec(tIdx >= 0) * 6;


%% ========================================================================
%  ALIGN DATA
%  ========================================================================

varKeys = fieldnames(varMap);

for iFile = 1:length(v)

    % Get idxPert
    idxPert = getFieldFromPath(v(iFile), pertMap);
    if isempty(idxPert), continue; end

    % Need reference nBins for this file to know where it ends
    nBinsRef = 0;
    for iVar = 1:length(varKeys)
        path = varMap.(varKeys{iVar});
        val = getFieldFromPath(v(iFile), path);

        if isnumeric(val) && ~isempty(val)
            nBinsRef = size(val, 2);
            break;
        elseif isstruct(val)
            flds = fieldnames(val);
            for f = 1:length(flds)
                dd = val.(flds{f});
                if isnumeric(dd) && ~isempty(dd) && size(dd, 2) > 1
                    nBinsRef = size(dd, 2);
                    break;
                end
            end
            if nBinsRef > 0, break; end
        end
    end

    if nBinsRef == 0, continue; end

    % Calculate extraction indices relative to THIS file's idxPert
    startIdx = idxPert - targetPre;
    endIdx = idxPert + targetPost - 1;

    % Helper function to pad and slice
    [padLeft, padRight, validStart, validEnd] = getPadding(startIdx, endIdx, nBinsRef);

    % --- Align All Variables ---
    for iKey = 1:length(varKeys)
        key = varKeys{iKey};
        path = varMap.(key);
        data = getFieldFromPath(v(iFile), path);

        if isempty(data), continue; end

        if isnumeric(data)
            % Check dimensions
            if size(data, 2) == nBinsRef
                dataAligned = sliceAndPad(data, validStart, validEnd, padLeft, padRight);
                v(iFile) = setFieldFromPath(v(iFile), path, dataAligned);
            end

        elseif isstruct(data)
            % Iterate fields
            flds = fieldnames(data);
            for iFld = 1:length(flds)
                fldName = flds{iFld};
                val = data.(fldName);

                % Check name first to prioritize time vector replacement
                if strcmp(fldName, 't') || strcmp(fldName, 'time')
                    data.(fldName) = tGlobal;

                elseif isnumeric(val) && ~isempty(val) && size(val, 2) == nBinsRef

                    valAligned = sliceAndPad(val, validStart, validEnd, padLeft, padRight);
                    data.(fldName) = valAligned;

                elseif isstruct(val) && isfield(val, 'idxPert')
                    % Update idxPert
                    data.(fldName).idxPert = globalIdxPert;
                    data.(fldName).tAlign = true;
                end
            end

            % Update struct in v
            v(iFile) = setFieldFromPath(v(iFile), path, data);
        end
    end
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [padLeft, padRight, validStart, validEnd] = getPadding(startIdx, endIdx, nBinsRef)
padLeft = 0;
padRight = 0;
validStart = startIdx;
validEnd = endIdx;

if startIdx < 1
    padLeft = 1 - startIdx;
    validStart = 1;
end

if endIdx > nBinsRef
    padRight = endIdx - nBinsRef;
    validEnd = nBinsRef;
end
end

function valAligned = sliceAndPad(data, validStart, validEnd, padLeft, padRight)
valSlice = data(:, validStart:validEnd);
nRows = size(data, 1);
valAligned = [nan(nRows, padLeft), valSlice, nan(nRows, padRight)];
end

function val = getFieldFromPath(s, path)
parts = strsplit(path, '.');
val = s;
for i = 1:length(parts)
    if isstruct(val) && isfield(val, parts{i})
        val = val.(parts{i});
    else
        val = [];
        return;
    end
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
    if ~isfield(s, parts{1})
        s.(parts{1}) = struct();
    end
    s.(parts{1}) = setFieldRecursive(s.(parts{1}), parts(2:end), val);
end
end
