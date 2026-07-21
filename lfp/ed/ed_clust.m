function [clustId, cInfo] = ed_clust(wv, tstamps, varargin)
% ED_CLUST Group candidate waveforms into shape types (PCA -> GMM).
%
%   [clustId, cInfo] = ED_CLUST(wv, tstamps, varargin)
%
%   SUMMARY:
%       Sorts a session's candidates by the shape of their LFP waveform, so
%       curation becomes a handful of decisions about TYPES instead of hundreds
%       about events.
%
%       It exists because a mean waveform lies. Average a set that holds
%       discharges, sharp waves and step artifacts and you get a curve that is
%       none of them, which is exactly how the discharges got buried the first
%       time this pipeline was calibrated. Split the set into shape clusters
%       first and each cluster's median is a real shape.
%
%       Deliberately NOT a classifier. It does not know which cluster is the
%       discharge - it only makes that question askable, and a human answers it
%       in ed_curate. Nothing here encodes the raMCU3/4/5 waveform, so a mouse
%       whose discharges look different still gets them in a cluster of their
%       own.
%
%       Preprocessing, every step measured on the 47 curated discharges
%       (dev/ed_clustSweep.m sweeps window x normalisation x features x count,
%       scored by how many events you must review to find 80% of them):
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
%           contrast: here the pool spans two orders of magnitude and the few
%           largest events take the principal components with them. Scale is
%           not lost - it comes back through the scalar measures below, ranked,
%           where it cannot dominate.
%         WINDOW +-50 ms, wider than their 10-50 ms. Also a different
%           contrast: a discharge and a sharp wave differ most in the DECAY (a
%           discharge is back to baseline in 50-100 ms, a sharp wave takes
%           300+), and a +-15 ms window cannot see it. Measured, +-50 and +-100
%           ms both beat +-15 and +-25 by a factor of three in review load.
%
%       Cluster count SCALES WITH THE POOL by default: k = 0.65*sqrt(n),
%       clamped. A fixed count cannot work across the range this sees - a pool
%       is 75 events under a strict filter and 8500 under a loose one, and 12
%       groups over 8500 leaves each one a mixture. Measured on raMCU3 (pool
%       8460, 10 curated discharges), the best cluster goes from 3% pure at
%       k=12 to 80% pure at k=60; the rule picks 60. On the old strict pools it
%       picks 11-14, which is where a fixed 12 had measured best.
%
%       Chosen over BIC, which maximises likelihood - not the objective - and
%       settled on ~7 everywhere.
%
%   INPUTS:
%       wv       - <mat>  [nEv x nSamp] per-event waveforms (edMaps.lfp).
%       tstamps  - <vec>  [1 x nSamp] window time base (s), from edMaps.
%       varargin - Parameter/Value:
%           'win'    - <vec> waveform window to cluster on (s). {[-0.05 0.05]}
%           'nPC'    - <num> principal components kept. {6}
%           'nClust' - <num> cluster count; empty scales it with the pool.
%                            {[]}
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
%                           largest cluster. NaN for an all-NaN row (an event
%                           too near a recording edge to have a waveform).
%       cInfo    - <struct> .nClust .score [nEv x nDim] .explained - what was
%                           fit, for the GUI and the record.
%
%   DEPENDENCIES:
%       pca, fitgmdist (Statistics and Machine Learning Toolbox).
%
%   HISTORY:
%       260721 created, replacing the threshold-knob curation. See
%              dev/ed_pipeline_rebuild.md.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'wv', @isnumeric);
addRequired(p, 'tstamps', @isnumeric);
addParameter(p, 'win', [-0.05 0.05], @isnumeric);
addParameter(p, 'nPC', 6, @isnumeric);
addParameter(p, 'nClust', [], @isnumeric);
addParameter(p, 'scalar', [], @isnumeric);
addParameter(p, 'detrend', 'edge', @ischar);
addParameter(p, 'norm', 'peak', @ischar);
addParameter(p, 'wSize', 0, @isnumeric);
parse(p, wv, tstamps, varargin{:});
win     = p.Results.win;
nPC     = p.Results.nPC;
nClust  = p.Results.nClust;
scalar  = p.Results.scalar;
flgDt   = lower(p.Results.detrend);
flgNorm = lower(p.Results.norm);
wSize   = p.Results.wSize;

nEv = size(wv, 1);
clustId = nan(nEv, 1);
cInfo = struct('nClust', 0, 'score', [], 'explained', []);

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

X = detrendWv(X, tstamps(iWin), flgDt);
[X, sz] = normWv(X, flgNorm);

% detrending costs two degrees of freedom, so X is rank-deficient by
% construction and pca says so on every call; the warning is expected, not a
% symptom
ws = warning('off', 'stats:pca:ColRankDefX');
oc = onCleanup(@() warning(ws));    % restored when the function exits
nPC = min(nPC, min(size(X)) - 1);
[~, score, ~, ~, explained] = pca(X, 'NumComponents', nPC);

% Scalar measures join the shape components on a comparable footing: rank
% first (fastZ and amp are heavy-tailed, and a raw one would set the metric by
% itself), then scale to the spread of the leading component.
if ~isempty(scalar)
    S = scalar(iOk, :);
    % A missing measure must not exile the event: .dur is NaN whenever no
    % half-amplitude crossing was found, which is ~16% of candidates, and
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

% SIZE, back as one explicit axis. Normalisation strips it from the waveform so
% the components can describe shape; wSize decides how much it then counts.
%
% On a LOG scale, because size is judged as a RATIO - a 3400 uV event against a
% 1000 uV one is the same kind of difference as 340 against 100, and the eye
% reads it that way while curating. The scalar block does carry .amp, but ranked
% (it is heavy-tailed), and a rank keeps only the ORDER: an event 3x larger than
% its neighbour and one 1.05x larger are the same distance apart. That is why
% raising amplitude's weight through the scalars could never work.
%
% Scaled to the spread of the leading shape component, so wSize = 1 means size
% counts for as much as the dominant shape axis and 0 restores shape-only.
if wSize > 0 && std(sz) > 0
    L = log10(max(sz, eps));
    L = (L - mean(L)) / max(std(L), eps);
    score = [score, L * std(score(:, 1)) * wSize];
end

% a GMM needs comfortably more events than dimensions; on a small pool keep
% only the leading features rather than failing
nDim = max(2, min(size(score, 2), floor(size(score, 1) / 5)));
score = score(:, 1 : nDim);

%% ========================================================================
%  CLUSTER
%  ========================================================================
% diagonal covariance and a regularisation floor: the clusters that matter are
% small (a discharge type can be a dozen events among hundreds), and a full
% covariance on a dozen points in 6 dimensions is singular.
gmOpt = {'CovarianceType', 'diagonal', 'RegularizationValue', 1e-6, ...
    'Replicates', 5, 'Options', statset('MaxIter', 500)};

% Fixed seed, restored on exit. fitgmdist starts from random centres, so
% without this the same pool and the same count give a different partition
% every call - and in the curation GUI that means pressing Re-cluster without
% changing anything reshuffles the groups the user just judged.
sRng = rng(0, 'twister');
ocRng = onCleanup(@() rng(sRng));

if isempty(nClust)
    nClust = round(0.65 * sqrt(size(score, 1)));
end
nClust = max(2, min(nClust, floor(size(score, 1) / 3)));
gm = fitgmdist(score, nClust, gmOpt{:});
lbl = cluster(gm, score);

% relabel by size, largest first, so a cluster index means something stable
% across a re-run and the GUI lists the bulk before the rare shapes
cnt = accumarray(lbl, 1, [nClust, 1]);
[~, ord] = sort(cnt, 'descend');
remap = zeros(nClust, 1);
remap(ord) = 1 : nClust;
clustId(iOk) = remap(lbl);

cInfo.nClust    = nClust;
cInfo.score     = nan(nEv, size(score, 2));
cInfo.score(iOk, :) = score;
cInfo.explained = explained(1 : min(nPC, numel(explained)));

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function X = detrendWv(X, t, met)
% Remove a per-event linear baseline, so an event riding on a slow deflection
% is described by its own shape.
%
%   'edge' fits the line on the FLANKS only - the samples beyond half the
%          window - so the event cannot tilt the baseline it is measured
%          against. This is the snipFromBinary convention.
%   'full' fits it over the whole window, the event included. Cheaper to say
%          and wrong in a specific way: a large asymmetric deflection drags
%          the line, and it drags it in a direction that depends on the
%          event's own polarity and asymmetry, so the residual tilt is
%          correlated with the shape being measured.
if strcmp(met, 'none'), return; end

t = t(:);
if strcmp(met, 'edge')
    iFit = abs(t) >= 0.5 * max(abs(t));
    if nnz(iFit) < 3, iFit = true(size(t)); end
else
    iFit = true(size(t));
end

t = t - mean(t(iFit));
A = [t(iFit), ones(nnz(iFit), 1)];
X = X - ([t, ones(numel(t), 1)] * (A \ X(:, iFit)'))';

end     % detrendWv


function [X, sz] = normWv(X, met)
% Put every event on a comparable scale, so the components describe SHAPE and
% not size. The pool spans two orders of magnitude, and unnormalised the few
% largest events take the principal components with them. Scale is not lost -
% it returns through the ranked scalar measures.
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
