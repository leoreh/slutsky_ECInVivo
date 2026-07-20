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
%       Detection and QA are separate concerns. Detection turns the signal into
%       candidate events (chMode..limDur, calibThr). QA turns events + per-event
%       metrics into an .accepted mask (.qa) - a pure "filters -> mask" spec that
%       ripp_gate applies. The same .qa is the headless default of ripp_curate
%       and the starting point its GUI seeds from; a user skips the GUI by
%       running with .qa as-is, or overrides it per mouse in the GUI.
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
%       .calibThr  - <log>  replace the fixed peak threshold with one calibrated
%                           to the 1/f noise floor (ripp_noiseFloor).
%       .targetFP  - <num>  target surrogate false-positive rate (events/s), used
%                           only when calibThr is true.
%       .qa        - <struct> the default QA filter (applied by ripp_gate):
%           .states - <vec>  vigilance-state indices to keep (AccuSleep order
%                            1=WAKE 2=QWAKE 3=LSLEEP 4=NREM 5=REM 6=N/REM
%                            7=BIN); [] = any state.
%           .ranges - <struct> per-event metric ranges [lo hi], one field per
%                            ripp per-event field. An event passes a metric if its
%                            value is inside the range or NaN (metric unavailable).
%                    The shipped default is deliberately PERMISSIVE - every state,
%                    no metric bounds - so detection keeps everything and the real
%                    gate is chosen per mouse in ripp_curate. Tighten it here only
%                    to change what an UNCURATED session falls back to. Note .emg
%                    is the robust z of evt_emgScore (median/MAD baseline), so a
%                    bound set here is on that scale.
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
%       260719b QA fields consolidated into .qa (the ripp_gate filter spec);
%               thrEmg/gainThr/nremOnly retired.
%       260720 help block corrected: it documented an active default filter
%              (.emg [-Inf 2], .spkGain [0 Inf], .states [2 3 4]) that the code
%              never set. The permissive default is the intended one - curation
%              owns the gate - so the doc was brought to the code, not vice versa.

if nargin < 1 || isempty(preset), preset = 'default'; end

% the shipping configuration (detection)
d = struct('name', 'default', 'chMode', 'tag', 'passband', [80 250], ...
    'detectMet', 3, 'zMet', 'nrem', ...
    'thr', [1 3.5 2 200 50], 'limDur', [15 300 20 10], ...
    'calibThr', false, 'targetFP', 0.05);

% the default QA filter (the ripp_gate spec): valid states, low EMG, MUA gate
d.qa = struct('states', [1 2 3 4 5 6 7]);
d.qa.ranges = struct('emg', [-Inf 1], 'spkGain', [1 Inf]);

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
