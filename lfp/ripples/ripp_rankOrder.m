function [rankMean, rankVar, timesFirst, timesLate] = ripp_rankOrder(spkTimes, rippTimes)
% RIPP_RANKORDER Calculates the normalized temporal rank of units within
% ripples.
%
%   INPUTS:
%       spkTimes      - (Cell) {N_units x 1} Spike times [s].
%       rippTimes     - (Mat)  [N_ripp x 2] Ripple start/end times [s].
%
%   OUTPUTS:
%       rankMean      - (Vec)  [N_units x 1] Mean normalized rank (0=Leader, 1=Follower).
%       rankVar       - (Vec)  [N_units x 1] Variance of rank order.
%       timesFirst    - (Cell) {N_units x 1} First spike in each ripple per unit.
%       timesLate     - (Cell) {N_units x 1} Subsequent spikes in each ripple per unit.
%
%   HISTORY:
%       Updated: 26 Jan 2026
%

% =========================================================================
%  CALCULATION
% =========================================================================

nUnits = length(spkTimes);

% Prepare Outputs
rankMean = nan(nUnits, 1);
rankVar  = nan(nUnits, 1);
timesFirst = cell(nUnits, 1);
timesLate  = cell(nUnits, 1);

% Flatten Spike Times for Vectorized Ops
allSpks  = [];
allUnits = [];
for iUnit = 1:nUnits
    if ~isempty(spkTimes{iUnit})
        allSpks  = [allSpks; spkTimes{iUnit}(:)]; %#ok<AGROW>
        allUnits = [allUnits; repmat(iUnit, length(spkTimes{iUnit}), 1)]; %#ok<AGROW>
    end
end

if isempty(allSpks)
    return;
end

% Map Spikes to Ripples
% Create edges for discretize: [Start1, End1, Start2, End2, ...]
rippEdges = reshape(rippTimes', [], 1);
binIdx    = discretize(allSpks, rippEdges);

% Keep only spikes inside ripple intervals (odd bins)
inRipp = mod(binIdx, 2) == 1;

relSpks  = allSpks(inRipp);
relUnits = allUnits(inRipp);
relRipp  = (binIdx(inRipp) + 1) / 2; % Convert bin index to Ripple ID

% Identify First Spike per Unit per Ripple
% Sort by RippleID then Timestamp
[~, sortIdx] = sortrows([relRipp, relSpks]);
srtdRipp  = relRipp(sortIdx);
srtdUnits = relUnits(sortIdx);
srtdSpks  = relSpks(sortIdx);

% Unique rows of [RippleID, UnitID] will return the first occurrence (lowest time)
[~, firstIdx] = unique([srtdRipp, srtdUnits], 'rows', 'first');

% Extract First Spikes
partRipp  = srtdRipp(firstIdx);
partUnits = srtdUnits(firstIdx);
partTimes = srtdSpks(firstIdx);

% Extract Late Spikes
lateMask = true(size(srtdSpks));
lateMask(firstIdx) = false;

lateUnits = srtdUnits(lateMask);
lateTimes = srtdSpks(lateMask);

% Pack Spike Times Outputs
% Re-accumulate into cells
for iSpk = 1:length(partTimes)
    uid = partUnits(iSpk);
    timesFirst{uid} = [timesFirst{uid}; partTimes(iSpk)];
end
for iSpk = 1:length(lateTimes)
    uid = lateUnits(iSpk);
    timesLate{uid} = [timesLate{uid}; lateTimes(iSpk)];
end

% Ensure sorted
timesFirst = cellfun(@sort, timesFirst, 'UniformOutput', false);
timesLate  = cellfun(@sort, timesLate,  'UniformOutput', false);

% =========================================================================
%  RANK CALCULATION
% =========================================================================

% Sort participants by RippleID (primary) and Time (secondary)
[~, rankSortIdx] = sortrows([partRipp, partTimes]);
partRipp  = partRipp(rankSortIdx);
partUnits = partUnits(rankSortIdx);

% Get number of participants per ripple
[~, ~, ic] = unique(partRipp);
countsPerRipp = accumarray(ic, 1);

% Calculate Rank (1-based index within each ripple group)
% Find start index of each ripple group in the sorted list
[~, grpStartIdx] = unique(partRipp, 'first');

% Expand group start index to every element
startIndices = grpStartIdx(ic);

% Rank = current_index - start_index + 1
ranks = (1:length(partRipp))' - startIndices + 1;

% Default to 0.5 (neutral) for single-participant ripples to avoid NaN
scores = ones(size(ranks)) * 0.5;

% Normalized Rank (0 to 1) -> (Rank - 1) / (Count - 1)
partCounts = countsPerRipp(ic);
maskMulti = partCounts > 1;
scores(maskMulti) = (ranks(maskMulti) - 1) ./ (partCounts(maskMulti) - 1);

% Aggregate per Unit
rankMean = accumarray(partUnits, scores, [nUnits 1], @mean, NaN);
rankVar  = accumarray(partUnits, scores, [nUnits 1], @var,  NaN);

end
