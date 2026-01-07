function frMat = cat_fr(basepaths)

%   CAT_FR concatenates firing rate matrices from different recording
%   sessions with respect to the absolute time of day. It stacks different
%   units from different sessions onto different rows of the final matrix.
%
%   Rows 1 : nUnits(1) are for the units from session 1,
%   Rows nUnits(1) + 1 : nUnits(1) + nUnits(2) are for session 2, etc.
%
% INPUT:
%   basepaths    (cell array) Full paths to recording folders.
%
% OUTPUT:
%   frMat        (matrix) Concatenated firing rate matrix [nTotalUnits x nTimeBins].

%% ========================================================================
%  INITIALIZE TIME INFO
%  ========================================================================

% Load 'fr' variable from all basepaths 
v = basepaths2vars('basepaths', basepaths, 'vars', {'fr'});

% Extract basenames and number of sessions
[~, basenames] = cellfun(@fileparts, basepaths, 'uni', false);
nPaths = length(basepaths);

% Define Reference Time (ZT0 = 09:00)
zt0 = guessDateTime('0900');

% Get recording start times
t_recStart = cellfun(@guessDateTime, basenames, 'uni', true);

% Define Experiment Start relative to ZT0 of the first day
t_expStart = t_recStart(1) - max([0, diff(timeofday([zt0, t_recStart(1)]))]);

% Get sampling interval (binsize)
binsize = v(1).fr.info.binsize;

% Calculate Start Index for the last session
recIdx_last = round(max([1, seconds(t_recStart(end) - t_expStart) / binsize]));

% Determine experiment end from the length of the fr of the last file
[~, nBins_last] = size(v(end).fr.strd);
expLen = recIdx_last + nBins_last - 1;

%% ========================================================================
%  LOAD AND CONCATENATE DATA
%  ========================================================================

frMat_chunks = {};

for iPath = 1 : nPaths

    currFr = v(iPath).fr.strd; % Assumes [nUnits x nBins]
    [nUnits, nBinsRec] = size(currFr);

    % Calculate Start Index for this session in the global matrix
    % Based on time difference from Experiment Start
    recIdx = round(max([1, seconds(t_recStart(iPath) - t_expStart) / binsize]));

    % Create padded matrix for this session's units
    sessionMat = nan(nUnits, expLen);

    % Determine valid range to fill
    validBins = min(nBinsRec, expLen - recIdx + 1);

    if validBins > 0
        % Place data into the global timeline
        sessionMat(:, recIdx : recIdx + validBins - 1) = currFr(:, 1 : validBins);
    end

    % Store chunk
    frMat_chunks{end+1} = sessionMat;
end

% Vertically concatenate all session chunks
if ~isempty(frMat_chunks)
    frMat = vertcat(frMat_chunks{:});
else
    frMat = [];
end

end

% EOF
