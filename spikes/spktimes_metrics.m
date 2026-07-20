function [st, stInfo] = spktimes_metrics(spktimes, bouts, varargin)
% SPKTIMES_METRICS Computes spike timing metrics inside one interval set.
%
%   [st, stInfo] = SPKTIMES_METRICS(spktimes, bouts, varargin)
%
%   SUMMARY:
%       Burstiness and irregularity metrics derived from the autocorrelogram
%       and the ISI distribution. Based in part on calc_ACG_metrics from cell
%       explorer. Two ACGs are computed: narrow (100 ms, 0.5 ms bins) and
%       wide (1 s, 1 ms bins).
%
%       The function takes ONE interval set and returns ONE value per unit.
%       Conditioning on states, days or epochs is done by spk_byCond, which
%       calls this repeatedly and labels the rows.
%
%       The intervals are treated as disjoint segments, not as a mask:
%
%       - ISIs are formed only between spikes of the same bout. Masking the
%         train and taking diff invents one interval per bout boundary.
%       - Bouts are pushed apart on a padded timeline before the ACG, so no
%         spike pair from different bouts can land inside the lag window.
%       - Each ACG lag is normalised by the number of spikes that had room
%         for a partner at that lag (nEff), not by the total spike count. A
%         spike near a bout edge cannot contribute at long lags; without this
%         a state made of short bouts shows a depressed acg baseline and an
%         inflated burst index purely from bout geometry.
%
%   INPUTS:
%       spktimes - (Cell)  {nUnits x 1} of spike times [s].
%       bouts    - (Mat)   [nBouts x 2] interval set [s]. Overlapping rows
%                          are consolidated. Pass [0 Inf] for the whole
%                          recording.
%       varargin - Parameter/Value pairs:
%           'fs'       - (Num)  Sampling frequency [Hz], for CCG. {20000}
%           'minSpks'  - (Num)  Metrics stay nan below this count. {100}
%           'basepath' - (Char) Recording path. {pwd}
%           'flgSave'  - (Log)  Save <basename>.st_metrics.mat. {false}.
%                               Off by default because the filename does not
%                               vary with bouts: saving from inside a
%                               spk_byCond map would leave only the last
%                               condition behind, wearing the whole-session
%                               name.
%
%   OUTPUTS:
%       st       - (Struct) Every field is [nUnits x ...] so the struct can
%                           go straight to struct2table:
%                    .nSpks     [nUnits x 1] spikes inside the bouts
%                    .dur       [nUnits x 1] total bout duration [s]
%                    .doublets  [nUnits x 1]
%                    .royer     [nUnits x 1]
%                    .royer2    [nUnits x 1]
%                    .lidor     [nUnits x 1]
%                    .mizuseki  [nUnits x 1]
%                    .cv        [nUnits x 1]
%                    .lv        [nUnits x 1]
%                    .acgNarrow [nUnits x nLagNarrow]
%                    .acgWide   [nUnits x nLagWide]
%       stInfo   - (Struct) Lag vectors [s] and parameters. Kept out of st so
%                           st stays table-able; folded back in on save.
%
%   DEPENDENCIES:
%       CCG, intervals.
%
%   HISTORY:
%       241121    LH
%       260719    one interval set in, one value out; segment-aware isi and
%                 acg; dropped the triple-exponential fit, lvr and cv2.
%
%   See also: SPK_BYCOND, BURST_STATS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addRequired(p, 'bouts', @isnumeric);
addParameter(p, 'fs', 20000, @isnumeric);
addParameter(p, 'minSpks', 100, @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, spktimes, bouts, varargin{:});
bouts    = p.Results.bouts;
fs       = p.Results.fs;
minSpks  = p.Results.minSpks;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% acg geometry. narrow resolves the refractory shoulder, wide reaches the
% 200-300 ms baseline that the royer index is normalised against
bnszNarrow = 0.0005;
durNarrow  = 0.1;
bnszWide   = 0.001;
durWide    = 1;

nUnits = numel(spktimes);

% an open bout is clipped to the last spike, so nEff below measures real
% recorded time rather than an infinite tail
if isempty(bouts), bouts = [0, Inf]; end
tMax = max(cellfun(@(x) max([0; x(:)]), spktimes));
bouts(isinf(bouts)) = tMax;

% consolidate so the bouts are sorted and disjoint; discretize below relies
% on one monotonic edge vector
bouts = intervals(bouts);
bouts = bouts.consolidate();
bouts = bouts.ints;

durBout = bouts(:, 2) - bouts(:, 1);
durTot  = sum(durBout);
edges   = reshape(bouts', [], 1);

% lag vectors, matching CCG's own geometry
halfNarrow = round(durNarrow / bnszNarrow / 2);
halfWide   = round(durWide / bnszWide / 2);
tNarrow    = (-halfNarrow : halfNarrow)' * bnszNarrow;
tWide      = (-halfWide : halfWide)' * bnszWide;

% initialize. nSpks and dur are always filled, so a unit that fails the
% count threshold reads as under-exposed rather than silently missing
st.nSpks     = zeros(nUnits, 1);
st.dur       = repmat(durTot, nUnits, 1);
st.doublets  = nan(nUnits, 1);
st.royer     = nan(nUnits, 1);
st.royer2    = nan(nUnits, 1);
st.lidor     = nan(nUnits, 1);
st.mizuseki  = nan(nUnits, 1);
st.cv        = nan(nUnits, 1);
st.lv        = nan(nUnits, 1);
st.acgNarrow = nan(nUnits, numel(tNarrow));
st.acgWide   = nan(nUnits, numel(tWide));


%% ========================================================================
%  RESTRICT AND CORRELATE
%  ========================================================================

% isi per unit, nan wherever the pair straddles a bout boundary
isiCell = cell(nUnits, 1);

for iUnit = 1 : nUnits

    % assign each spike to a bout. an odd bin index means inside a bout, an
    % even one means the gap between two, nan means outside the set
    spkT = sort(spktimes{iUnit}(:));
    iBin = discretize(spkT, edges);
    inBout = ~isnan(iBin) & mod(iBin, 2) == 1;

    spkT    = spkT(inBout);
    boutIdx = (iBin(inBout) + 1) / 2;
    nSpks   = numel(spkT);

    st.nSpks(iUnit) = nSpks;
    if nSpks < minSpks
        continue
    end

    % isis, with cross-bout pairs marked rather than removed so that lv can
    % still tell which isis were consecutive in the original train
    isiAll = diff(spkT);
    isiAll(boutIdx(1 : end - 1) ~= boutIdx(2 : end)) = NaN;
    isiCell{iUnit} = isiAll;

    % push the bouts apart so the gap between any two exceeds the widest lag
    % window, then correlate the whole train in one call
    pad = durWide;
    offset = [0; cumsum(durBout(1 : end - 1) + pad)];
    tPad = spkT - bouts(boutIdx, 1) + offset(boutIdx);

    acgN = CCG(tPad, ones(nSpks, 1), 'binSize', bnszNarrow, ...
        'duration', durNarrow, 'norm', 'counts', 'Fs', 1 / fs);
    acgW = CCG(tPad, ones(nSpks, 1), 'binSize', bnszWide, ...
        'duration', durWide, 'norm', 'counts', 'Fs', 1 / fs);

    % rate normalisation, with the edge correction
    dLeft  = spkT - bouts(boutIdx, 1);
    dRight = bouts(boutIdx, 2) - spkT;

    st.acgNarrow(iUnit, :) = acgN(:)' ./ ...
        (nEffLag(dLeft, dRight, tNarrow) * bnszNarrow);
    st.acgWide(iUnit, :) = acgW(:)' ./ ...
        (nEffLag(dLeft, dRight, tWide) * bnszWide);
end


%% ========================================================================
%  METRICS
%  ========================================================================

% bin offsets from the zero-lag bin. integers rather than lag comparisons so
% the windows cannot drift on a float rounding
iZeroN = halfNarrow + 1;
iZeroW = halfWide + 1;

% royer divides by the 200-300 ms baseline, which can be zero for a sparse
% unit. one pseudocount, half the smallest positive baseline in the
% population, keeps the ratio finite without inventing structure
bslAll = mean(st.acgWide(:, iZeroW + 200 : iZeroW + 300), 2, 'omitnan');
pseudoCnt = min(bslAll(bslAll > 0));
if isempty(pseudoCnt)
    pseudoCnt = 1;
end
pseudoCnt = pseudoCnt / 2;

for iUnit = 1 : nUnits

    if st.nSpks(iUnit) < minSpks
        continue
    end

    acgN   = st.acgNarrow(iUnit, :);
    acgW   = st.acgWide(iUnit, :);
    isiAll = isiCell{iUnit};

    % burstiness ---------------------------------------------------------

    % doublets: peak of the 2.5-8 ms bins over the mean of the 8-11.5 ms bins
    st.doublets(iUnit) = max(acgN(iZeroN + 5 : iZeroN + 16)) / ...
        mean(acgN(iZeroN + 16 : iZeroN + 23));

    % royer 2012: mean of the 3-5 ms bins over the mean of 200-300 ms
    brstMean = mean(acgW(iZeroW + 3 : iZeroW + 5));
    bslMean  = mean(acgW(iZeroW + 200 : iZeroW + 300));
    if bslMean <= 0
        st.royer(iUnit) = (brstMean + pseudoCnt) / (bslMean + pseudoCnt);
    else
        st.royer(iUnit) = brstMean / bslMean;
    end

    % royer2: peak of 0-10 ms against the 40-50 ms baseline, normalised by
    % whichever of the two is larger so the index stays bounded
    brstPeak = max(acgW(iZeroW : iZeroW + 10));
    bslPeak  = mean(acgW(iZeroW + 40 : iZeroW + 50));
    if brstPeak > bslPeak
        st.royer2(iUnit) = (brstPeak - bslPeak) / brstPeak;
    else
        st.royer2(iUnit) = (brstPeak - bslPeak) / bslPeak;
    end

    % lidor: 2-10 ms against 35-50 ms as a contrast index
    brstSum = sum(acgN(iZeroN + 5 : iZeroN + 19));
    bslSum  = sum(acgN(iZeroN + 71 : iZeroN + 99));
    st.lidor(iUnit) = (brstSum - bslSum) / (brstSum + bslSum);

    % mizuseki 2011: fraction of spikes flanked by an isi below 6 ms on
    % either side. a cross-bout isi is inf, so a bout boundary can neither
    % create a burst spike nor hide one
    isiNbr = [inf; isiAll; inf];
    isiNbr(isnan(isiNbr)) = inf;
    st.mizuseki(iUnit) = mean(isiNbr(1 : end - 1) < 0.006 | ...
        isiNbr(2 : end) < 0.006);

    % firing irregularity ------------------------------------------------

    % cv: shinomoto 2003. the cross-bout isis are gone, so this is now the
    % within-state dispersion rather than the gap structure between bouts
    isi = isiAll(~isnan(isiAll));
    if numel(isi) < 2
        continue
    end
    st.cv(iUnit) = std(isi) / mean(isi);

    % lv: shinomoto 2003, kobayashi 2019. local, so a pair drops out unless
    % both isis were consecutive and inside the same bout
    isi1 = isiAll(1 : end - 1);
    isi2 = isiAll(2 : end);
    lvTerm = 3 * (isi1 - isi2).^2 ./ (isi1 + isi2).^2;
    st.lv(iUnit) = mean(lvTerm, 'omitnan');
end


%% ========================================================================
%  FINALIZE
%  ========================================================================

stInfo.bouts      = bouts;
stInfo.dur        = durTot;
stInfo.tNarrow    = tNarrow;
stInfo.tWide      = tWide;
stInfo.bnszNarrow = bnszNarrow;
stInfo.bnszWide   = bnszWide;
stInfo.minSpks    = minSpks;
stInfo.fs         = fs;
stInfo.runtime    = datetime("now");

if flgSave
    [~, basename] = fileparts(basepath);
    stFile = fullfile(basepath, [basename, '.st_metrics.mat']);
    backup_file(stFile);

    % one variable per file, named st, so that basepaths2vars keeps
    % resolving 'st_metrics' to v(i).st
    sVar.st = st;
    sVar.st.info = stInfo;
    save(stFile, '-struct', 'sVar')
end

end     % MAIN


%% ========================================================================
%  LOCALS
%  ========================================================================

function nEff = nEffLag(dLeft, dRight, tLag)
% Number of reference spikes with room for a partner at each lag. At lag tau
% only a spike at least |tau| from the relevant bout edge could have had one;
% counting every spike instead tapers the correlogram toward long lags in
% proportion to how short the bouts are.

tLag = tLag(:)';
nEff = nan(1, numel(tLag));

% lags run ascending and symmetric about zero. the negative side is mirrored
% so that both threshold vectors ascend, as histcounts edges must
iPos = tLag >= 0;
nEff(iPos)  = cntAtLeast(dRight, tLag(iPos));
nEff(~iPos) = flip(cntAtLeast(dLeft, flip(-tLag(~iPos))));

% a lag no spike can reach carries no information; nan keeps it out of the
% metrics instead of dividing by zero
nEff(nEff == 0) = NaN;

end     % nEffLag


function n = cntAtLeast(d, thr)
% Count elements of d that are >= each threshold. thr must ascend. Sentinel
% edges bracket the data so that histcounts sees a finite, strictly
% increasing edge vector and every element lands in exactly one bin.

lo = min([d(:); thr(:)]) - 1;
hi = max([d(:); thr(:)]) + 1;

n = numel(d) - cumsum(histcounts(d, [lo, thr, hi]));
n = n(1 : end - 1);

end     % EOF
