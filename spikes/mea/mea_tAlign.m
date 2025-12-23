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
    % Access reference variable using helper logic (simple field access)
    % Support nested struct access if needed (e.g. 'mea.fr')
    try
        dataRef = getFieldFromPath(v(iFile), refPath);
    catch
        continue; % Skip if ref var missing
    end

    if isempty(dataRef), continue; end

    % Check consistency of binSize
    if isempty(binSize)
        binSize = dataRef.info.binSize;
    elseif dataRef.info.binSize ~= binSize
        warning('Bin size mismatch in file %d. Expected %f, got %f.', iFile, binSize, dataRef.info.binSize);
    end

    idxPert = dataRef.info.idxPert;
    nBins = size(dataRef.fr, 2); % Assumes .fr field exists in reference

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

    % Calculates padding required for this file
    % We need to access the reference again to know this file's specific offset
    try
        dataRef = getFieldFromPath(v(iFile), refPath);
    catch
        continue;
    end

    if isempty(dataRef), continue; end

    idxPert = dataRef.info.idxPert;
    nBinsRef = size(dataRef.fr, 2);

    nPre = idxPert - 1;
    padPre = maxPre - nPre;

    currentPost = nBinsRef - idxPert;
    padPost = maxPost - currentPost;


    % Now iterate over all variables to align
    for iKey = 1:length(varKeys)
        key = varKeys{iKey};
        path = varMap.(key);

        try
            data = getFieldFromPath(v(iFile), path);
        catch
            continue;
        end

        if isempty(data), continue; end

        % Iterate fields of the variable struct
        flds = fieldnames(data);
        for iFld = 1:length(flds)
            fldName = flds{iFld};
            val = data.(fldName);

            % Check name first to prioritize time vector replacement
            if strcmp(fldName, 't') || strcmp(fldName, 'time')
                % Replace time vector
                data.(fldName) = tGlobal;

                % Check if it matches reference bins (numeric matrix)
                % Dims: (nUnits x nBins) or (1 x nBins)
            elseif isnumeric(val) && ~isempty(val) && size(val, 2) == nBinsRef

                nRows = size(val, 1);

                % Pad
                valAligned = [nan(nRows, padPre), val, nan(nRows, padPost)];
                data.(fldName) = valAligned;

            elseif strcmp(fldName, 'info')
                % Update info if it's the reference variable
                if strcmp(key, refVar)
                    data.info.idxPert = globalIdxPert;
                    data.info.tAlign = true;
                end
            end
        end

        % Put data back into v
        v(iFile) = setFieldFromPath(v(iFile), path, data);

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
