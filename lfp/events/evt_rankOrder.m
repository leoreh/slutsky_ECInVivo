function [rankMean, rankVar, timesFirst, timesLate] = evt_rankOrder(spkTimes, evtTimes)
% EVT_RANKORDER Calculates the normalized temporal rank of units within
% events.
%
%   INPUTS:
%       spkTimes      - (Cell) {N_units x 1} Spike times [s].
%       evtTimes     - (Mat)  [N_ripp x 2] Event start/end times [s].
%
%   OUTPUTS:
%       rankMean      - (Vec)  [N_units x 1] Mean normalized rank (0=Leader, 1=Follower).
%       rankVar       - (Vec)  [N_units x 1] Variance of rank order.
%       timesFirst    - (Cell) {N_units x 1} First spike in each event per unit.
%       timesLate     - (Cell) {N_units x 1} Subsequent spikes in each event per unit.
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

% Map Spikes to Events
% Create edges for discretize: [Start1, End1, Start2, End2, ...]
evtEdges = reshape(evtTimes', [], 1);
binIdx    = discretize(allSpks, evtEdges);

% Keep only spikes inside event intervals (odd bins)
inEvt = mod(binIdx, 2) == 1;

relSpks  = allSpks(inEvt);
relUnits = allUnits(inEvt);
relEvt  = (binIdx(inEvt) + 1) / 2; % Convert bin index to Event ID

% Identify First Spike per Unit per Event
% Sort by EventID then Timestamp
[~, sortIdx] = sortrows([relEvt, relSpks]);
srtdEvt  = relEvt(sortIdx);
srtdUnits = relUnits(sortIdx);
srtdSpks  = relSpks(sortIdx);

% Unique rows of [EventID, UnitID] will return the first occurrence (lowest time)
[~, firstIdx] = unique([srtdEvt, srtdUnits], 'rows', 'first');

% Extract First Spikes
partEvt  = srtdEvt(firstIdx);
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

% Sort participants by EventID (primary) and Time (secondary)
[~, rankSortIdx] = sortrows([partEvt, partTimes]);
partEvt  = partEvt(rankSortIdx);
partUnits = partUnits(rankSortIdx);

% Get number of participants per event
[~, ~, ic] = unique(partEvt);
countsPerRipp = accumarray(ic, 1);

% Calculate Rank (1-based index within each event group)
% Find start index of each event group in the sorted list
[~, grpStartIdx] = unique(partEvt, 'first');

% Expand group start index to every element
startIndices = grpStartIdx(ic);

% Rank = current_index - start_index + 1
ranks = (1:length(partEvt))' - startIndices + 1;

% Default to 0.5 (neutral) for single-participant events to avoid NaN
scores = ones(size(ranks)) * 0.5;

% Normalized Rank (0 to 1) -> (Rank - 1) / (Count - 1)
partCounts = countsPerRipp(ic);
maskMulti = partCounts > 1;
scores(maskMulti) = (ranks(maskMulti) - 1) ./ (partCounts(maskMulti) - 1);

% Aggregate per Unit
rankMean = accumarray(partUnits, scores, [nUnits 1], @mean, NaN);
rankVar  = accumarray(partUnits, scores, [nUnits 1], @var,  NaN);

end
