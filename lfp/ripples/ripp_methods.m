function met = ripp_methods(preset)
% RIPP_METHODS Detection + QA configuration (met) for the ripple pipeline.
%
%   met = RIPP_METHODS(preset)
%
%   SUMMARY:
%       The one place that defines a ripple "method": a struct carrying every
%       detection and quality-assurance choice, so ripp_wrapper (and ripp_screen)
%       take one met instead of a long argument list. 'default' returns the
%       single shipping configuration; 'screen' returns an array of methods to
%       compare (the shipping detector vs the same detector with its threshold
%       calibrated to the 1/f noise floor). Edit a field here, not in the wrapper.
%
%   FIELDS:
%       .name      - <char> short id, used as the report column / provenance.
%       .chMode    - <char> 'tag' = follow ripp.info.rippCh (else best channel);
%                           'best' = force a fresh best-NREM-channel pick.
%       .passband  - <vec>  band-pass [lo hi] (Hz).
%       .detectMet - <num>  ripp_sigPrep detection signal (3 = TEO).
%       .zMet      - <char> ripp_sigPrep threshold reference ('nrem').
%       .thr       - <vec>  ripp_times thresholds [start peak cont max minCont].
%       .limDur    - <vec>  ripp_times duration limits [min max inter minCont] (ms).
%       .thrEmg    - <num>  QA: reject an event whose EMG z exceeds this.
%       .gainThr   - <num>  QA: keep an event only if its MUA gain exceeds this
%                           (the false-positive gate; tune interactively with
%                           ripp_gateGui; 0 keeps any positive-gain event).
%       .nremOnly  - <log>  QA: accept only NREM events (else valid states 2/3/4).
%       .calibThr  - <log>  replace the fixed peak threshold with one calibrated
%                           to the 1/f noise floor (ripp_noiseFloor).
%       .targetFP  - <num>  target surrogate false-positive rate (events/s), used
%                           only when calibThr is true.
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

if nargin < 1 || isempty(preset), preset = 'default'; end

% the shipping configuration
d = struct('name', 'default', 'chMode', 'tag', 'passband', [80 250], ...
    'detectMet', 3, 'zMet', 'nrem', ...
    'thr', [1 3.5 2 200 50], 'limDur', [15 300 20 10], ...
    'thrEmg', 2, 'gainThr', 0, 'nremOnly', false, ...
    'calibThr', false, 'targetFP', 0.05);

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
