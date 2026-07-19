function m = ripp_screenMethods()
% RIPP_SCREENMETHODS Detection methods compared by ripp_screen.
%
%   m = RIPP_SCREENMETHODS()
%
%   SUMMARY:
%       The one place to define, edit, add, or drop a method. Set the shared
%       defaults once, then give each method a block that overrides only what
%       differs. To add a method copy a block, bump the index, change a field;
%       to drop one delete its block. Methods that share chMode, passband,
%       detectMet, and zMet share one filtered signal in ripp_screen (it is
%       prepared once), so a pair that differs only in the threshold costs
%       almost nothing extra.
%
%   FIELDS (per method):
%       .name      - <char> short id, used as the report column and table key.
%       .chMode    - <char> 'tag'  = average the ripp channel (evt_rippCh, the
%                           shipping pipeline); 'best' = single channel with the
%                           most ripple-band power in NREM.
%       .passband  - <vec>  band-pass [lo hi] (Hz).
%       .detectMet - <num>  ripp_sigPrep detection signal (3 = TEO).
%       .zMet      - <char> ripp_sigPrep threshold reference ('nrem').
%       .thr       - <vec>  ripp_times thresholds [start peak cont max minCont].
%       .limDur    - <vec>  ripp_times duration limits [min max inter minCont] (ms).
%       .calibThr  - <log>  replace the fixed peak threshold with one calibrated
%                           to the 1/f noise floor (ripp_noiseFloor).
%       .targetFP  - <num>  target surrogate false-positive rate (events/s),
%                           used only when calibThr is true.
%
%   HISTORY:
%       260716 detection-review parameter screen.
%       260717 recast as current vs 1/f-calibrated (fooof) threshold.

% shared defaults - the shipping pipeline; a block overrides only what differs
d = struct('name', '', 'chMode', 'tag', 'passband', [80 250], ...
    'detectMet', 3, 'zMet', 'nrem', ...
    'thr', [1 3.5 2 200 50], 'limDur', [15 300 20 10], ...
    'calibThr', false, 'targetFP', 0.05);

% current - fixed SD threshold (the shipping pipeline). The baseline column.
m(1)      = d;
m(1).name = 'current';

% fooof - same signal, threshold calibrated to the recording's 1/f noise floor
m(2)          = d;
m(2).name     = 'fooof';
m(2).calibThr = true;

end     % EOF
