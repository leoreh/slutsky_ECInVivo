function [lfp, emg, fs] = ripp_sigLoad(basepath, varargin)
% RIPP_SIGLOAD Load the ripple detection channel + EMG, windowed.
%
%   [lfp, emg, fs] = RIPP_SIGLOAD(basepath, varargin)
%
%   SUMMARY:
%       Loads the signals the ripple pipeline detects on, mirroring ed_sigLoad.
%       Reads the ripple LFP channel from <basename>.lfp (averaged if several,
%       in microvolts) over the analysis window, and the matching EMG trace from
%       <basename>.sleep_sig.mat. The EMG is windowed to the same span and, as a
%       fail-safe, dropped when its length does not match the LFP (a sampling
%       mismatch between the two files) - which relaxes the EMG QA criterion
%       rather than scoring events against a misaligned trace.
%
%   INPUTS:
%       basepath - (Char) Session directory.
%       varargin - Parameter/Value pairs:
%           'win'      - (Vec)    Window [start end] (s). {[0 Inf]}
%           'session'  - (Struct) Session metadata (loaded if empty).
%           'basename' - (Char)   File stem. {folder name}
%           'rippCh'   - (Num)    Detection channel. {ripp.info.rippCh or best}
%           'bit2uv'   - (Num)    Conversion factor. {auto}
%
%   OUTPUTS:
%       lfp      - (Vec) [n x 1] Windowed ripple channel (microvolts).
%       emg      - (Vec) [n x 1] Windowed EMG, or [] on mismatch / absence.
%       fs       - (Num) LFP sampling frequency [Hz].
%
%   DEPENDENCIES:
%       ripp_pickCh, evt_loadCh, basepaths2vars (only when session is empty).
%
%   HISTORY:
%       Created: 260706 (parity with ed_sigLoad; via the shared loader).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'bit2uv', [], @isnumeric);
parse(p, basepath, varargin{:});

basepath = p.Results.basepath;
win      = p.Results.win;
session  = p.Results.session;
rippCh   = p.Results.rippCh;
bit2uv   = p.Results.bit2uv;

basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end

% Session metadata (needed for channel, rates, conversion)
if isempty(session)
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'});
    session = v.session;
end

%% ========================================================================
%  RIPPLE CHANNEL
%  ========================================================================
if isempty(rippCh)
    rippCh = ripp_pickCh(basepath, 'basename', basename, 'session', session, ...
        'win', win);
end
[lfp, fs] = evt_loadCh(basepath, basename, session, rippCh, win, bit2uv);

%% ========================================================================
%  EMG  (windowed to the LFP; dropped on length mismatch)
%  ========================================================================
emg = [];
sigFile = fullfile(basepath, [basename, '.sleep_sig.mat']);
if isfile(sigFile)
    S = load(sigFile, 'emg');
    if isfield(S, 'emg'), emg = S.emg(:); end
end

if ~isempty(emg)
    s1 = round(win(1) * fs) + 1;
    if isinf(win(2))
        emg = emg(s1:end);
    else
        s2 = min(length(emg), round(win(2) * fs));
        emg = emg(s1:s2);
    end
    if length(emg) ~= length(lfp)
        warning('ripp_sigLoad:emgFit', ...
            'EMG (%d) and LFP (%d) lengths differ; skipping EMG QA.', ...
            length(emg), length(lfp));
        emg = [];
    end
end

end     % EOF
