function rippSpks = ripp_spks(spkTimes, rippTimes, ctrlTimes, varargin)
% RIPP_SPKS Analyzes spiking rate modulation during ripples.
%
%   rippSpks = RIPP_SPKS(spkTimes, rippTimes, ctrlTimes, varargin)
%
%   SUMMARY:
%       Calculates scalar modulation metrics comparing Ripple vs Control periods.
%       1. Instantaneous Firing Rates (FR) per event.
%       2. Mean FR differences (RippleMu vs ControlMu).
%       3. Z-Scored Gain.
%       4. Statistical Significance (Wilcoxon Sign-Rank / Rank-Sum).
%       5. Rank Order (Mean and Variance).
%
%   INPUTS:
%       spkTimes    - (Cell) {N_units x 1} Spike times [s].
%       rippTimes   - (Mat)  [N x 2] Ripple start/end times [s].
%       ctrlTimes   - (Mat)  [N x 2] Control start/end times [s].
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgSave'  - (Log)  Save output? (Default: true).
%
%   OUTPUTS:
%       rippSpks    - (Struct) Stats structure with [N_units x 1] fields:
%           .frRipp    - Mean Firing Rate during Ripples (Hz).
%           .frRand    - Mean Firing Rate during Control (Hz).
%           .frZ       - Z-scored modulation ((Ripp - Ctrl) / StdCtrl).
%           .frMod     - Modulation Index ((Ripp - Ctrl) / (Ripp + Ctrl)).
%           .pVal      - P-value from significance test.
%           .h0        - Boolean significance flag (p < 0.05).
%           .pFire     - Probability of participation (fraction of ripples with spikes).
%           .spkCount  - Mean number of spikes per ripple.
%           .rankMean  - Mean normalized rank order (0=Leader, 1=Follower).
%           .rankVar   - Variance of rank order.
%           .frActive  - Mean firing rate during active ripples (Hypersynchrony).
%           .fracUnits - (nRipples x 1) Fraction of units participating in each ripple.
%
%   DEPENDENCIES:
%       times2rate.
%
%   HISTORY:
%       Updated: 23 Jan 2026

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
parse(p, spkTimes, rippTimes, ctrlTimes, varargin{:});

basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippSpks.mat']);

% =========================================================================
%  CALCULATION
% =========================================================================

nUnits = length(spkTimes);

% Calculate Durations
durRipp = (rippTimes(:, 2) - rippTimes(:, 1))'; % [1 x nRipp]
durCtrl = (ctrlTimes(:, 2) - ctrlTimes(:, 1))'; % [1 x nRand]

% Get Spike COUNTS
rippCounts = times2rate(spkTimes, 'winCalc', rippTimes, 'binsize', Inf, 'c2r', false);
ctrlCounts = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf, 'c2r', false);

% Calculate RATES (Hz)
rippRates = rippCounts ./ durRipp;
ctrlRates = ctrlCounts ./ durCtrl;

% METRICS
% Mean Rates
frRipp = mean(rippRates, 2, 'omitnan');
frRand = mean(ctrlRates, 2, 'omitnan');
sdRand = std(ctrlRates, [], 2, 'omitnan');

% Spike Count Stats
spkCount = mean(rippCounts, 2, 'omitnan');
pFire    = mean(rippCounts > 0, 2, 'omitnan');

% Modulation Metrics
frZ = (frRipp - frRand) ./ sdRand;
frZ(sdRand < eps) = NaN; % Avoid infs

frSum = frRipp + frRand;
frMod = (frRipp - frRand) ./ frSum;
frMod(frSum < eps) = NaN;

% Conditional Rate
% Captures "Hypersynchrony" independent of "Reliability" (Intensity when
% Active). Essentially equal to frRipp ./ pFire.
% Prepare matrix with NaNs where rate is 0 (inactive)
activeMat = rippRates;
activeMat(activeMat == 0) = NaN;
frActive = mean(activeMat, 2, 'omitnan');

% =========================================================================
%  STATISTICAL SIGNIFICANCE
% =========================================================================

pVal = nan(nUnits, 1);
flgPair = size(rippRates, 2) == size(ctrlRates, 2);

% Pre-calculate nanflag to speed up loop skip
goodIdx = ~all(isnan(rippRates), 2) & sum(rippCounts, 2) > 0;

for iUnit = 1:nUnits
    if ~goodIdx(iUnit), continue; end

    if flgPair
        pVal(iUnit) = signrank(rippRates(iUnit, :), ctrlRates(iUnit, :));
    else
        pVal(iUnit) = ranksum(rippRates(iUnit, :), ctrlRates(iUnit, :));
    end
end
h0 = pVal < 0.05;

% =========================================================================
%  RANK ORDER
% =========================================================================

% Flatten Spike Times for Vectorized Ops
allSpks  = [];
allUnits = [];
for iUnit = 1:nUnits
    if ~isempty(spkTimes{iUnit})
        allSpks  = [allSpks; spkTimes{iUnit}(:)]; %#ok<AGROW>
        allUnits = [allUnits; repmat(iUnit, length(spkTimes{iUnit}), 1)]; %#ok<AGROW>
    end
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
% Sort by RippleID then Timestamp to ensure 'unique' picks the first one
[~, sortIdx] = sortrows([relRipp, relSpks]);
srtdRipp  = relRipp(sortIdx);
srtdUnits = relUnits(sortIdx);
srtdSpks  = relSpks(sortIdx);

% Unique rows of [RippleID, UnitID] will return the first occurrence (lowest time)
% because we just sorted by time.
[~, firstIdx] = unique([srtdRipp, srtdUnits], 'rows', 'first');

% Keep only first spikes
partRipp  = srtdRipp(firstIdx);
partUnits = srtdUnits(firstIdx);
partTimes = srtdSpks(firstIdx);

% Calculate Ranks
% Sort again by RippleID (primary) and Time (secondary) to establish rank
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

% Normalize Rank (0 to 1) -> Rank / Count
partCounts = countsPerRipp(ic);
scores(partCounts > 1) = (ranks(partCounts > 1) - 1) ./ (partCounts(partCounts > 1) - 1);

% Aggregate per Unit
rankMean = accumarray(partUnits, scores, [nUnits 1], @mean, NaN);
rankVar  = accumarray(partUnits, scores, [nUnits 1], @var,  NaN);

% Fraction of Units Participating per Ripple
nRipp = size(rippTimes, 1);
countsAllRipp = accumarray(partRipp, 1, [nRipp 1]);
fracUnits = countsAllRipp / nUnits;

% =========================================================================
%  OUTPUT & SAVE
% =========================================================================

% Pack results
rippSpks.frRipp    = frRipp;
rippSpks.frRand    = frRand;
rippSpks.frZ       = frZ;
rippSpks.frMod     = frMod;
rippSpks.spkCount  = spkCount;
rippSpks.pFire     = pFire;
rippSpks.pVal      = pVal;
rippSpks.h0    = h0;
rippSpks.rankMean  = rankMean;
rippSpks.rankVar   = rankVar;
rippSpks.frActive  = frActive;
rippSpks.fracUnits = fracUnits;


if flgSave
    save(savefile, 'rippSpks', '-v7.3');
end

end