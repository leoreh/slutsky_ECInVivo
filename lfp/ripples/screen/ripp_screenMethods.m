function m = ripp_screenMethods()
% RIPP_SCREENMETHODS Detection methods compared by ripp_screen.
%
%   m = RIPP_SCREENMETHODS()
%
%   The one place to define, edit, add, or drop a method. Set the shared
%   defaults once, then give each method its own block and override only what
%   differs. To add a method: copy a block, bump the index, change a field.
%   To drop one: delete its block. Nothing else in the screen changes.
%
%   FIELDS (per method):
%       .name      - <char> short id, used as the report column.
%       .chMode    - <char> 'tag'  = average the channelTags.Ripple channels
%                           (the shipping pipeline), 'best' = the single
%                           channel with the most ripple-band power in NREM.
%       .passband  - <vec>  band-pass [lo hi] (Hz).
%       .detectMet - <num>  ripp_sigPrep detection signal (3 = TEO).
%       .zMet      - <char> ripp_sigPrep threshold reference:
%                           'nrem'   = z over NREM, 'nremBg' = z over NREM with
%                           candidate ripples removed (signal-independent).
%       .thr       - <vec>  ripp_times thresholds [start peak cont max minCont].
%       .limDur    - <vec>  ripp_times duration limits [min max inter minCont] (ms).
%
%   HISTORY:
%       260716 detection-review parameter screen.

% shared defaults - a method's block overrides only what differs
d = struct('name', '', 'chMode', 'best', 'passband', [120 220], ...
    'detectMet', 3, 'zMet', 'nrem', ...
    'thr', [1 3.5 2 200 50], 'limDur', [15 300 20 10]);

% current - the shipping pipeline (tagged channels, wide band). The baseline.
m(1)          = d;
m(1).name     = 'current';
m(1).chMode   = 'tag';
m(1).passband = [80 250];

% narrow - the proposed fix: best channel, ripple band only.
m(2)          = d;
m(2).name     = 'narrow';

% narrowAbs - narrow, but threshold off a signal-independent baseline, to see
% whether per-animal z hides a group difference.
m(3)          = d;
m(3).name     = 'narrowAbs';
m(3).zMet     = 'nremBg';

end     % EOF
