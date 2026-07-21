function met = ripp_methods(preset)
% RIPP_METHODS Detection + QA configuration (met) for the ripple pipeline.
%
%   met = RIPP_METHODS(preset)
%
%   SUMMARY:
%       The one place that defines a ripple "method": a struct carrying every
%       detection choice plus the default QA filter, so the pipeline stages take
%       one met instead of a long argument list. 'default' returns the single
%       shipping configuration; 'screen' returns an array of methods to compare
%       (the shipping detector vs the same detector with its threshold calibrated
%       to the 1/f noise floor). Edit a field here, not in the stages.
%
%       Three concerns, kept apart:
%         DETECTION turns the signal into candidate events (chMode..limDur,
%           calibThr). Permissive on purpose - a candidate costs a row.
%         .qa is a NOISE FILTER, not a definition. It asks only the two
%           questions no ripple can fail: was the muscle quiet, did units fire.
%           Its job is to hand the clustering a workable pool, and to be the
%           mask an uncurated session falls back to. evt_gate applies it.
%         .clust groups the survivors by waveform shape, and a human accepts
%           whole shapes in ripp_curate. That is where "is this a ripple" is
%           decided, which is why no shape criterion lives in .qa.
%
%   FIELDS:
%       .name      - <char> short id, used as the report column / provenance.
%       .chMode    - <char> 'tag' = follow ripp.info.rippCh (else best channel);
%                           'best' = force a fresh best-NREM-channel pick.
%       .passband  - <vec>  band-pass [lo hi] (Hz).
%       .detectMet - <num>  ripp_sigPrep detection signal (3 = TEO).
%       .zMet      - <char> ripp_sigPrep threshold reference ('nrem').
%       .otlThr    - <num>  gross-artifact threshold in robust SDs (lfp_artifacts).
%                           Flagged samples are dropped from the 'nrem' baseline
%                           only - never from detection. One constant for every
%                           mouse; it self-scales, landing ~1-5 mV. 8 sits on a
%                           plateau (8-12): below 6 the mask reaches real events,
%                           above 15 it stops catching movement steps.
%       .thr       - <vec>  ripp_times thresholds [start peak cont max minCont].
%       .limDur    - <vec>  ripp_times duration limits [min max inter minCont] (ms).
%       .calibThr  - <log>  replace the fixed peak threshold with one calibrated
%                           to the 1/f noise floor (ripp_noiseFloor).
%       .targetFP  - <num>  target surrogate false-positive rate (events/s), used
%                           only when calibThr is true.
%       .qa        - <struct> the NOISE FILTER that sizes the pool handed to
%                            the clustering, and the mask an UNCURATED session
%                            falls back to (applied by evt_gate):
%           .states - <vec>  vigilance-state indices to keep (AccuSleep order
%                            1=WAKE 2=QWAKE 3=LSLEEP 4=NREM 5=REM 6=N/REM
%                            7=BIN); [] = any state. In the GUI this is a SCOPE,
%                            not a verdict - the state boxes start from this
%                            list and can be flipped back and forth.
%           .ranges - <struct> per-event metric ranges [lo hi], one field per
%                            ripp per-event field. An event passes a metric if
%                            its value is inside the range or NaN (metric
%                            unavailable). Each FINITE bound becomes one knob in
%                            the curation GUI, so a metric added here appears
%                            there with no other change. Note .emg is the robust
%                            z of evt_emgScore (median/MAD baseline), so a bound
%                            set here is on that scale.
%                    These two ask only what no ripple can fail. What a ripple
%                    LOOKS like is not asked here - that is .clust's business,
%                    answered per mouse in ripp_curate.
%       .clust     - <struct> evt_clust arguments (.win .nPC .nClust .nFit
%                            .detrend .norm .wSize), plus .scalar (the per-event
%                            fields joining the shape components) and .nView
%                            (the per-cluster row cap for the curation view).
%
%   INPUTS:
%       preset - <char> 'default' (one met) | 'screen' (comparison array).
%                       {default}
%
%   OUTPUT:
%       met    - <struct> one configuration, or an array of them.
%
%   HISTORY:
%       260716 detection-review parameter screen (as ripp_screenMethods).
%       260719 folded into the pipeline as the met config; 'current' -> 'default'.
%       260719b QA fields consolidated into .qa (the evt_gate filter spec);
%               thrEmg/gainThr/nremOnly retired.
%       260720 help block corrected: it documented an active default filter
%              (.emg [-Inf 2], .spkGain [0 Inf], .states [2 3 4]) that the code
%              never set. The permissive default is the intended one - curation
%              owns the gate - so the doc was brought to the code, not vice versa.
%       260722 .clust added: curation became waveform clustering (ripp_curate ->
%              evt_curate), so .qa stopped being the answer and became the noise
%              filter that sizes the pool handed to it.

if nargin < 1 || isempty(preset), preset = 'default'; end

% the shipping configuration (detection)
d = struct('name', 'default', 'chMode', 'tag', 'passband', [80 250], ...
    'detectMet', 3, 'zMet', 'nrem', 'otlThr', 8, ...
    'thr', [1 3.5 2 200 50], 'limDur', [15 300 20 10], ...
    'calibThr', false, 'targetFP', 0.05);

% The noise filter that sizes the POOL (the evt_gate spec). Two criteria, each
% blind to the other's failure mode and neither describing a shape:
%   emg      quiet muscle -> not a movement transient
%   spkGain  units fired  -> not a bystander deflection
% Every state is kept; the GUI's state boxes start from this list.
d.qa = struct('states', [1 2 3 4 5 6 7]);
d.qa.ranges = struct('emg', [-Inf 1], 'spkGain', [1 Inf]);

% Waveform clustering (evt_clust), the stage that decides what a ripple is.
%
% WIN is +-30 ms, narrower than the ED pipeline's +-50: a ripple separates from
% a step or a spike-bleed transient in the oscillation itself, while a discharge
% separates from a sharp wave in the DECAY, which needs the wider window.
%
% NCLUST is an explicit 20 rather than the sqrt rule evt_clust falls back to.
% The rule is built for finding a RARE shape (10 discharges among 8460
% candidates, where 12 groups leave every group a mixture); a ripple pool's
% contaminant is a whole population, so a coarse partition separates it, and
% 0.65*sqrt(30000) would ask a human to read 110 tiles. Set the GUI box to 0 for
% the rule when a cluster is obviously mixed.
%
% NFIT and NVIEW bound the cost at ripple scale (pools of 8k-35k against an ED
% pool of 3k-8k): the mixture is estimated on 8000 events and the rest are
% projected and assigned, and the view is handed at most 300 rows per cluster.
% Neither touches the mask - every event is labelled, and the counts on the
% checkboxes are the true ones.
d.clust = struct('win', [-0.03 0.03], 'nPC', 6, 'nClust', 20, ...
    'detrend', 'edge', 'norm', 'peak', 'wSize', 1, ...
    'nFit', 8000, 'nView', 300);
d.clust.scalar = {'peakProm', 'freqPeak', 'amp', 'dur', 'emg', 'spkGain'};

switch preset
    case 'default'
        met = d;

    case 'screen'
        % the shipping detector vs the same signal with a 1/f-calibrated
        % threshold (the only difference is calibThr)
        met(1)          = d;
        met(2)          = d;
        met(2).name     = 'fooof';
        met(2).calibThr = true;

    otherwise
        error('ripp_methods:preset', 'unknown preset "%s"', preset);
end

end     % EOF
