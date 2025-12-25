function [v, tGlobal] = mea_tAlign(v, varMap, pertMap, flgEdge)
% MEA_TALIGN Aligns firing rate matrices and other variables to perturbation onset.
%
%   [v, t] = MEA_TALIGN(v, varMap, pertMap, flgEdge) receives a struct array
%   'v' and aligns specified variables such that the perturbation index
%   (idxPert) is identical across all files. It also returns a global time
%   vector 't'.
%
%   INPUTS:
%       v           - (struct array) Data structure containing variables to align.
%       varMap      - (struct) Mapping of variable names to align.
%                     Fields are arbitrary names, values are paths in 'v'.
%                     Example: varMap.fr = 'mea.fr'; varMap.dyn = 'dyn';
%       pertMap     - (char) Path to the perturbation index variable in 'v'.
%                     Example: 'fr.info.idxPert'.
%       flgEdge     - (logical, optional)
%                     If true (default): Intersection alignment. Clips data to
%                     the latest start and earliest end times across all files.
%                     Values before the latest start or after the earliest end
%                     are discarded.
%                     If false: Union alignment. Pads shorter files with NaNs
%                     to match the earliest start and latest end times.
%
%   OUTPUTS:
%       v           - (struct array) Updated struct with aligned fields.
%       tGlobal     - (vector) Global time vector [seconds].
%
%   See also: MEA_FRPREP, BASEPATHS2VARS

if nargin < 4 || isempty(flgEdge)
    flgEdge = true;
end

%% ========================================================================
%  ANALYZE ALIGNMENT DIMENSIONS
%  ========================================================================

maxPre = 0;
maxPost = 0;
minPre = inf;
minPost = inf;
binSize = [];

% Determine the pre- and post-perturbation durations for all files
for iFile = 1:length(v)

    % Access perturbation index
    idxPert = getFieldFromPath(v(iFile), pertMap);
    if isempty(idxPert)
        warning('Could not find idxPert for file %d at %s', iFile, pertMap);
        continue;
    end

    % Attempt to find binSize (assume sibling of idxPert)
    if isempty(binSize)
        parts = strsplit(pertMap, '.');
        if length(parts) > 1
            infoPath = strjoin(parts(1:end-1), '.');
            iInfo = getFieldFromPath(v(iFile), infoPath);
            if isstruct(iInfo) && isfield(iInfo, 'binSize')
                binSize = iInfo.binSize;
            end
        end
    end

    % Determine nBins from the first valid numeric variable in varMap
    nBins = 0;
    varKeys = fieldnames(varMap);
    for k = 1:length(varKeys)
        path = varMap.(varKeys{k});
        val = getFieldFromPath(v(iFile), path);

        if isnumeric(val) && ~isempty(val)
            nBins = size(val, 2);
            break;
        elseif isstruct(val)
            % Try to find existing numeric fields inside struct (e.g. .fr)
            flds = fieldnames(val);
            for f=1:length(flds)
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
    % pre-pert bins: 1 to (idxPert - 1) -> Count is idxPert - 1
    % post-pert bins: idxPert to nBins  -> Count is nBins - idxPert + 1 (Wait)
    % Convention: idxPert is the index of the perturbation bin (t=0 or first post bin).
    % Usually idxPert is the first bin of the "post" epoch.
    % So Pre bins = idxPert - 1.
    % Post bins = nBins - idxPert + 1. (Including the perturbation bin itself).
    % Let's stick to simple "bins before" and "bins including/after".

    currentPre = idxPert - 1;
    currentPost = nBins - idxPert + 1; % Include the pert bin in post count for convenience

    if currentPre > maxPre, maxPre = currentPre; end
    if currentPost > maxPost, maxPost = currentPost; end

    if currentPre < minPre, minPre = currentPre; end
    if currentPost < minPost, minPost = currentPost; end
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

if isempty(binSize)
    tGlobal = [];
    warning('No valid binSize found. Cannot create global time vector.');
else
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
end


%% ========================================================================
%  ALIGN DATA
%  ========================================================================

varKeys = fieldnames(varMap);

for iFile = 1:length(v)

    % Get idxPert
    idxPert = getFieldFromPath(v(iFile), pertMap);
    if isempty(idxPert), continue; end

    % Need reference nBins for this file to know where it ends
    % (Re-detect nBinsRef same way as above)
    nBinsRef = 0;
    for k = 1:length(varKeys)
        path = varMap.(varKeys{k});
        val = getFieldFromPath(v(iFile), path);
        if isnumeric(val) && ~isempty(val)
            nBinsRef = size(val, 2);
            break;
        elseif isstruct(val)
            flds = fieldnames(val);
            for f=1:length(flds)
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
    % We want [idxPert - targetPre, ..., idxPert + targetPost - 1]
    % (Note: targetPost included bin 0 so -1 if we just sum counts)

    startIdx = idxPert - targetPre;
    endIdx = idxPert + targetPost - 1;

    % If flgEdge=false (UNION), these indices might be Out Of Bounds.
    % We need to pad.

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


    % --- Align All Variables ---
    for iKey = 1:length(varKeys)
        key = varKeys{iKey};
        path = varMap.(key);
        data = getFieldFromPath(v(iFile), path);

        if isempty(data), continue; end

        if isnumeric(data)
            % If numeric, check dimensions
            if size(data, 2) == nBinsRef
                % Slice valid part
                valSlice = data(:, validStart:validEnd);

                % Pad
                nRows = size(data, 1);
                valAligned = [nan(nRows, padLeft), valSlice, nan(nRows, padRight)];

                v(iFile) = setFieldFromPath(v(iFile), path, valAligned);
            end

        elseif isstruct(data)
            % If struct, iterate fields
            flds = fieldnames(data);
            for iFld = 1:length(flds)
                fldName = flds{iFld};
                val = data.(fldName);

                % Check name first to prioritize time vector replacement
                if strcmp(fldName, 't') || strcmp(fldName, 'time')
                    data.(fldName) = tGlobal;

                elseif isnumeric(val) && ~isempty(val) && size(val, 2) == nBinsRef

                    valSlice = val(:, validStart:validEnd);
                    nRows = size(val, 1);
                    valAligned = [nan(nRows, padLeft), valSlice, nan(nRows, padRight)];

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

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

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
