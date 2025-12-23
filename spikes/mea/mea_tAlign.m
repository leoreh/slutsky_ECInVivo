function v = mea_tAlign(basepaths, varargin)
% MEA_TALIGN Aligns firing rate matrices from multiple recordings to perturbation onset.
%
%   v = MEA_TALIGN(BASEPATHS) loads 'fr' structures from the specified
%   paths, aligns them such that the perturbation index (idxPert) is
%   identical across all recordings (by NaN-padding), and updates the time
%   vector to a global timeframe.
%
%   INPUTS:
%       basepaths   - (cell) Cell array of directory paths containing .fr.mat files.
%
%   OPTIONAL (Key-Value Pairs):
%       flgSave     - (logical) Whether to save the updated struct back to disk. {false}
%
%   OUTPUTS:
%       v           - (struct array) aligned 'fr' structures.
%
%   See also: MEA_FRPREP, BASEPATHS2VARS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'basepaths', @iscell);
addParameter(p, 'flgSave', false, @islogical);

parse(p, basepaths, varargin{:});
flgSave = p.Results.flgSave;


%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load fr structures using basepaths2vars
vStruct = basepaths2vars('basepaths', basepaths, 'vars', "fr");

% Extract the 'fr' field from the loaded structure to work with
% We keep vStruct to preserve path information for saving
if isempty(vStruct)
    v = [];
    return;
end


%% ========================================================================
%  ANALYZE ALIGNMENT
%  ========================================================================

maxPre = 0;
maxPost = 0;
binSize = [];

% Determine the maximum pre- and post-perturbation durations
for iFile = 1:length(vStruct)
    fr = vStruct(iFile).fr;
    if isempty(fr), continue; end

    % Check consistency of binSize
    if isempty(binSize)
        binSize = fr.info.binSize;
    elseif fr.info.binSize ~= binSize
        warning('Bin size mismatch in file %d. Expected %f, got %f.', iFile, binSize, fr.info.binSize);
    end

    idxPert = fr.info.idxPert;
    nBins = size(fr.fr, 2);

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

v = []; % Container for output structs

for iFile = 1:length(vStruct)
    fr = vStruct(iFile).fr;

    if isempty(fr)
        frAligned = [];
    else
        % Calculate padding required
        idxPert = fr.info.idxPert;
        nBins = size(fr.fr, 2);

        nPre = idxPert - 1;
        % padPre is amount to add before current start
        padPre = maxPre - nPre;

        % padPost is amount to add after current end to reach maxPost
        % Current post is (nBins - idxPert)
        % We want maxPost
        currentPost = nBins - idxPert;
        padPost = maxPost - currentPost;

        nUnits = size(fr.fr, 1);

        % Pad matrices with NaNs
        fr.fr = [nan(nUnits, padPre), fr.fr, nan(nUnits, padPost)];
        fr.frOrig = [nan(nUnits, padPre), fr.frOrig, nan(nUnits, padPost)];

        % Update info and time
        fr.info.idxPert = globalIdxPert;
        fr.t = tGlobal;
        fr.info.tAlign = true;

        frAligned = fr;

        % Save if requested
        if flgSave
            basepath = basepaths{iFile};
            [~, basename] = fileparts(basepath);
            saveFile = fullfile(basepath, [basename, '.fr.mat']);

            % We need to save the variable 'fr'
            save(saveFile, 'fr', '-v7.3');
        end
    end

    % Store in output
    if iFile == 1
        v = frAligned;
    else
        v(iFile) = frAligned;
    end
end

end     % EOF
