function [clustId, cInfo] = evt_clust(wv, tstamps, varargin)
% EVT_CLUST Group event waveforms into shape types (PCA -> GMM).
%
%   [clustId, cInfo] = EVT_CLUST(wv, tstamps, varargin)
%
%   SUMMARY:
%       Sorts a session's events by the shape of their LFP waveform, so
%       curation becomes a handful of decisions about TYPES instead of
%       thousands about events. Shared by both event pipelines (evt_curate).
%
%       It exists because a mean waveform lies. Average a set that holds
%       discharges, sharp waves and step artifacts and you get a curve that is
%       none of them, which is exactly how the discharges got buried the first
%       time the ED pipeline was calibrated. Split the set into shape clusters
%       first and each cluster's median is a real shape.
%
%       Deliberately NOT a classifier. It does not know which cluster is the
%       event of interest - it only makes that question askable, and a human
%       answers it in evt_curate. Nothing here encodes a particular waveform,
%       so a mouse whose events look different still gets them in a cluster of
%       their own.
%
%       Preprocessing, every step measured on the 47 curated discharges
%       (lfp/ed/dev/ed_clustSweep.m sweeps window x normalisation x features x
%       count, scored by how many events you must review to find 80% of them):
%
%         DETREND, fitting the line on the FLANKS, which removes the slow
%           deflection the event happens to sit on without letting the event
%           tilt its own baseline. Fitting over the whole window instead leaves
%           a residual tilt that depends on the event's polarity and asymmetry
%           - a distortion correlated with the very shape being measured. The
%           47 curated discharges cannot separate the two (dev/ed_alignSweep.m:
%           the ranking flips with how the sweep is aggregated), so this one is
%           chosen on the argument, not the score. It is also the convention
%           snipFromBinary has used for spike waveforms since 2020.
%         NORMALISE each waveform to unit peak. Maslarova et al. 2025 keep
%           absolute amplitude for ripple-versus-IED, but that is a different
%           contrast: an ED pool spans two orders of magnitude and the few
%           largest events take the principal components with them. Scale is
%           not lost - it comes back through .wSize and the scalar measures.
%         WINDOW is the caller's (met.clust.win), because the shapes being told
%           apart differ per pipeline. A discharge separates from a sharp wave
%           in the DECAY, which needs +-50 ms; a ripple separates from a step
%           artifact in the oscillation itself, which is over in +-30 ms.
%
%       Cluster count SCALES WITH THE POOL when nClust is empty:
%       k = 0.65*sqrt(n), clamped. A fixed count cannot work across the range
%       an ED pool sees - 75 events under a strict filter and 8500 under a
%       loose one, and 12 groups over 8500 leaves each one a mixture. Measured
%       on raMCU3 (pool 8460, 10 curated discharges), the best cluster goes
%       from 3% pure at k=12 to 80% pure at k=60; the rule picks 60. Chosen
%       over BIC, which maximises likelihood - not the objective - and settled
%       on ~7 everywhere.
%
%       A RIPPLE pool is a different problem and passes an explicit count. Its
%       contaminant is a whole population rather than a rare shape, so a coarse
%       partition suffices, and the sqrt rule over 30k events would ask a human
%       to read 110 tiles.
%
%       NFIT bounds the cost. Both the PCA basis and the mixture are estimates
%       of the pool's shape, and an estimate does not improve once the sample
%       is large: fitting on nFit events and then PROJECTING and ASSIGNING all
%       of them gives the same partition for a fraction of the work. A 30k-event
%       ripple pool is minutes of fitgmdist at full size and seconds at 8000.
%       Empty = fit on everything, which is the ED default and leaves that
%       pipeline's partition bit-identical.
%
%   INPUTS:
%       wv       - <mat>  [nEv x nSamp] per-event waveforms (evtMaps.lfp).
%       tstamps  - <vec>  [1 x nSamp] window time base (s), from evt_maps.
%       varargin - Parameter/Value:
%           'win'    - <vec> waveform window to cluster on (s). {[-0.05 0.05]}
%           'nPC'    - <num> principal components kept. {6}
%           'nClust' - <num> cluster count; empty scales it with the pool.
%                            {[]}
%           'nFit'   - <num> events the PCA + GMM are fit on (drawn at random,
%                            deterministically); empty fits on all. {[]}
%           'detrend'- <char> per-event baseline removal: 'edge' fits the line
%                            on the flanks only, 'full' over the whole window
%                            (the event included, so it tilts its own
%                            baseline), 'none'. {'edge'}
%           'norm'   - <char> per-event scaling: 'peak' (unit peak, L-inf),
%                            'l2', 'none'. {'peak'}
%           'wSize'  - <num>  weight of log10(size) as an extra axis, in units
%                            of the leading shape component's spread. 0 is
%                            shape only. Effectively on/off - a diagonal GMM
%                            re-fits variance per dimension, so the value
%                            barely changes the partition. {0}
%           'scalar' - <mat> [nEv x nFeat] extra per-event measures to cluster
%                            on alongside the shape components. Each is
%                            rank-normalised, so a heavy-tailed one cannot
%                            dominate the geometry. Passing them measured
%                            better than shape alone on every mouse and every
%                            window, so the pipeline always does. {[]}
%
%   OUTPUTS:
%       clustId  - <vec>    [nEv x 1] cluster index, ordered so 1 is the
%                           TALLEST cluster (largest median peak height). NaN
%                           for an all-NaN row (an event too near a recording
%                           edge to have a waveform).
%       cInfo    - <struct> .nClust .score [nEv x nDim] .explained .nFit -
%                           what was fit, for the GUI and the record.
%
%   DEPENDENCIES:
%       evt_detrend; pca, fitgmdist (Statistics and Machine Learning Toolbox).
%
%   HISTORY:
%       260721 created as lfp/ed/ed_clust, replacing threshold-knob curation.
%              See lfp/ed/dev/ed_pipeline_rebuild.md.
%       260722 moved to lfp/events as evt_clust and shared with the ripple
%              pipeline (the body was already event-agnostic; the same move
%              ripp_gate -> evt_gate made). Gained 'nFit', because a ripple
%              pool is 5-20x an ED pool and fitgmdist is superlinear in it.
%       260722b clusters ordered by peak height, not population size, so the
%              curation GUI can lay them on an amplitude continuum.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'wv', @isnumeric);
addRequired(p, 'tstamps', @isnumeric);
addParameter(p, 'win', [-0.05 0.05], @isnumeric);
addParameter(p, 'nPC', 6, @isnumeric);
addParameter(p, 'nClust', [], @isnumeric);
addParameter(p, 'nFit', [], @isnumeric);
addParameter(p, 'scalar', [], @isnumeric);
addParameter(p, 'detrend', 'edge', @ischar);
addParameter(p, 'norm', 'peak', @ischar);
addParameter(p, 'wSize', 0, @isnumeric);
parse(p, wv, tstamps, varargin{:});
win     = p.Results.win;
nPC     = p.Results.nPC;
nClust  = p.Results.nClust;
nFit    = p.Results.nFit;
scalar  = p.Results.scalar;
flgDt   = lower(p.Results.detrend);
flgNorm = lower(p.Results.norm);
wSize   = p.Results.wSize;

nEv = size(wv, 1);
clustId = nan(nEv, 1);
cInfo = struct('nClust', 0, 'score', [], 'explained', [], 'nFit', 0);

%% ========================================================================
%  FEATURES
%  ========================================================================
iWin = tstamps >= win(1) & tstamps <= win(2);
X = double(wv(:, iWin));

% an event at a recording edge has a NaN waveform and cannot be clustered;
% it keeps a NaN label and the caller treats that as its own group
iOk = all(isfinite(X), 2);
% below this a "type" is not a meaningful thing to fit, and the caller treats
% an all-NaN labelling as "too few to cluster"
MINEV = 20;
X = X(iOk, :);
if size(X, 1) < MINEV || size(X, 2) < 3
    return;
end

X = evt_detrend(X, tstamps(iWin), flgDt);
pk = max(abs(X), [], 2);            % peak height, for ordering the clusters
[X, sz] = normWv(X, flgNorm);

% Fixed seed, restored on exit. Both the fit subsample and fitgmdist's random
% starts draw from it, so the same pool and count give the same partition every
% call - in the curation GUI, pressing Re-cluster without changing anything must
% be a no-op rather than a reshuffle of the groups just judged. With nFit empty
% nothing is drawn here, so the stream fitgmdist sees is untouched.
sRng = rng(0, 'twister');
ocRng = onCleanup(@() rng(sRng));

% The rows the estimates are fit on. Everything else is projected and assigned,
% which is O(n) rather than O(n * k * iter).
% A cap below MINEV is read as NO cap, not as "fit on 20 events": a caller
% asking for a fit smaller than the floor a fit needs has misconfigured it, and
% silently estimating a 20-cluster mixture from 20 events would return a
% partition that looks fine and means nothing.
nOk = size(X, 1);
iFit = (1 : nOk)';
if ~isempty(nFit) && isfinite(nFit) && nFit >= MINEV && nFit < nOk
    iFit = sort(randperm(nOk, round(nFit)))';
end
nF = numel(iFit);

% detrending costs two degrees of freedom, so X is rank-deficient by
% construction and pca says so on every call; the warning is expected, not a
% symptom
ws = warning('off', 'stats:pca:ColRankDefX');
oc = onCleanup(@() warning(ws));    % restored when the function exits
nPC = min(nPC, min(nF, size(X, 2)) - 1);
[coeff, ~, ~, ~, explained, mu] = pca(X(iFit, :), 'NumComponents', nPC);
score = (X - mu) * coeff;           % identical to pca's own score when nF = nOk

% SIZE, back as one explicit axis. Normalisation strips it from the waveform so
% the components can describe shape; wSize decides how much it then counts.
%
% On a LOG scale, because size is judged as a RATIO - a 3400 uV event against a
% 1000 uV one is the same kind of difference as 340 against 100, and the eye
% reads it that way while curating. The scalar block below does carry .amp, but
% ranked (it is heavy-tailed), and a rank keeps only the ORDER: an event 3x
% larger than its neighbour and one 1.05x larger are the same distance apart.
% That is why raising amplitude's weight through the scalars could never work.
%
% Scaled to the spread of the leading shape component, so wSize = 1 means size
% counts for as much as the dominant shape axis and 0 restores shape-only.
%
% It sits BEFORE the scalars because a scarce-dimension pool truncates from the
% right (see nDim): shape first, then size, then the derived measures.
if wSize > 0 && std(sz) > 0
    L = log10(max(sz, eps));
    L = (L - mean(L)) / max(std(L), eps);
    score = [score, L * std(score(:, 1)) * wSize];
end

% Scalar measures join the shape components on a comparable footing: rank
% first (fastZ and amp are heavy-tailed, and a raw one would set the metric by
% itself), then scale to the spread of the leading component. Ranked over ALL
% events, not the fit subsample, so a projected event lands on the same axis as
% a fitted one.
if ~isempty(scalar)
    S = scalar(iOk, :);
    % A missing measure must not exile the event: .dur is NaN whenever no
    % half-amplitude crossing was found, which is ~16% of ED candidates, and
    % dropping those would make them unclusterable and therefore permanently
    % unacceptable in the GUI. The waveform is what drives the grouping, so an
    % absent scalar takes its column median - neutral in the ranking.
    for iCol = 1 : size(S, 2)
        bad = ~isfinite(S(:, iCol));
        if ~any(bad), continue; end
        fill = median(S(~bad, iCol));
        if ~isfinite(fill), fill = 0; end    % the whole column is missing
        S(bad, iCol) = fill;
    end
    S = (tiedrank(S) - 0.5) ./ size(S, 1);
    score = [score, (S - 0.5) * std(score(:, 1)) * 2];
end

% a GMM needs comfortably more events than dimensions; on a small pool keep
% only the leading features rather than failing. Sized on the FIT set, which is
% what the mixture actually sees.
nDim = max(2, min(size(score, 2), floor(nF / 5)));
score = score(:, 1 : nDim);

%% ========================================================================
%  CLUSTER
%  ========================================================================
% diagonal covariance and a regularisation floor: the clusters that matter are
% small (a discharge type can be a dozen events among hundreds), and a full
% covariance on a dozen points in 6 dimensions is singular.
gmOpt = {'CovarianceType', 'diagonal', 'RegularizationValue', 1e-6, ...
    'Replicates', 5, 'Options', statset('MaxIter', 500)};

if isempty(nClust)
    nClust = round(0.65 * sqrt(nOk));
end
nClust = max(2, min(nClust, floor(nF / 3)));
gm = fitgmdist(score(iFit, :), nClust, gmOpt{:});
lbl = cluster(gm, score);

% relabel by peak height, tallest first, so the GUI lists and tiles the
% clusters on a scannable amplitude continuum instead of by arbitrary
% population size - the tall step artifacts a ripple pool carries then sit
% together at the top, next to the tall genuine ripples. The order is a display
% convenience only; a saved partition keeps its own labels regardless.
h = accumarray(lbl, pk, [nClust, 1], @median);
[~, ord] = sort(h, 'descend');
remap = zeros(nClust, 1);
remap(ord) = 1 : nClust;
clustId(iOk) = remap(lbl);

cInfo.nClust    = nClust;
cInfo.score     = nan(nEv, size(score, 2));
cInfo.score(iOk, :) = score;
cInfo.explained = explained(1 : min(nPC, numel(explained)));
cInfo.nFit      = nF;

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function [X, sz] = normWv(X, met)
% Put every event on a comparable scale, so the components describe SHAPE and
% not size. An ED pool spans two orders of magnitude, and unnormalised the few
% largest events take the principal components with them. Scale is not lost -
% it returns through .wSize and the ranked scalar measures.
%
%   'peak' divides by max|x| (L-inf). One noisy sample sets the whole scale.
%   'l2'   divides by the norm over the window, which no single sample can
%          dominate - but which a long tail inflates, so a wider window makes
%          a slow event look smaller.
switch met
    case 'none'
        sz = ones(size(X, 1), 1);
        return
    case 'l2'
        sz = vecnorm(X, 2, 2);
    otherwise
        sz = max(abs(X), [], 2);
end
sz(sz == 0) = 1;
X = X ./ sz;

end     % normWv
