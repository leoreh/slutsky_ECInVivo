function [sig, emg, fs, edCh] = ed_sigLoad(basepath, varargin)
% ED_SIGLOAD Load the ED detection signal + EMG, windowed.
%
%   [sig, emg, fs, edCh] = ED_SIGLOAD(basepath, varargin)
%
%   SUMMARY:
%       Loads the signals the ED pipeline detects on, mirroring ripp_sigLoad.
%       Two channel modes, chosen by met.chMode:
%       - 'eeg'  (default) the assembled sleep signal <basename>.sleep_sig.mat
%                 eeg field. This is the only source that exists on every
%                 session that was sleep-scored, including recordings kept
%                 without a session.mat or a binary .lfp (the EA cohort), so it
%                 is what the pipeline falls back to and what the legacy IED
%                 detection ran on.
%       - 'ripp'  the raw channel ripple detection ran on, read from
%                 <basename>.lfp via the shared evt_loadCh. Use it to put EDs
%                 and ripples on the same electrode in one mouse; it needs a
%                 session.mat and a binary .lfp.
%
%       The EMG always comes from sleep_sig and is windowed to the same span.
%       As a fail-safe it is dropped when its length does not match the
%       detection signal - a rate mismatch between the .lfp and sleep_sig -
%       which leaves the EMG metric NaN (and therefore permissive in evt_gate)
%       rather than scoring events against a misaligned trace.
%
%   INPUTS:
%       basepath - (Char) Session directory.
%       varargin - Parameter/Value pairs:
%           'basename' - (Char)   File stem. {folder name}
%           'chMode'   - (Char)   'eeg' | 'ripp'. {'eeg'}
%           'edCh'     - (Num)    Explicit channel for 'ripp'; else resolved by
%                                 ripp_pickCh. Ignored by 'eeg'. {[]}
%           'win'      - (Vec)    Window [start end] (s). {[0 Inf]}
%           'session'  - (Struct) Session metadata (loaded if empty & needed).
%           'bit2uv'   - (Num)    Conversion for 'ripp'. {auto}
%
%   OUTPUTS:
%       sig      - (Vec)  [n x 1] Windowed detection signal.
%       emg      - (Vec)  [n x 1] Windowed EMG, or [] on mismatch / absence.
%       fs       - (Num)  Sampling frequency of SIG [Hz].
%       edCh     - (Num)  Channel actually used ([] for the 'eeg' mode).
%
%   DEPENDENCIES:
%       evt_loadCh, ripp_pickCh, basepaths2vars (only for the 'ripp' mode).
%
%   HISTORY:
%       Created: 260622
%       Updated: 260706 (lfp branch via the shared evt_loadCh).
%       Updated: 260720 (rebuilt for the staged pipeline: 'sigSource' became
%                met.chMode with 'eeg' the default on both sides - the wrapper
%                and the loader used to disagree - and the four unread outputs
%                (emgRms, specAdapter, sSig) were dropped with the old ed_gui
%                that consumed them. The EMG is now dropped on a length
%                mismatch instead of being scored at the wrong rate.)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'chMode', 'eeg', @(x) any(strcmpi(x, {'eeg', 'ripp'})));
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'bit2uv', [], @isnumeric);
parse(p, basepath, varargin{:});

basepath = p.Results.basepath;
chMode   = lower(p.Results.chMode);
edCh     = p.Results.edCh;
win      = p.Results.win;
session  = p.Results.session;
bit2uv   = p.Results.bit2uv;

basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
sigFile = fullfile(basepath, [basename, '.sleep_sig.mat']);

%% ========================================================================
%  DETECTION SIGNAL
%  ========================================================================
switch chMode
    case 'eeg'
        if ~isfile(sigFile)
            error('ed_sigLoad:noFile', 'sleep_sig.mat not found: %s', sigFile);
        end
        S   = load(sigFile, 'eeg', 'fs');
        sig = S.eeg(:);
        fs  = S.fs;
        edCh = [];

    case 'ripp'
        if isempty(session)
            v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'});
            session = v.session;
        end
        if isempty(edCh)
            edCh = ripp_pickCh(basepath, 'basename', basename, ...
                'session', session, 'win', win);
        end
        [sig, fs] = evt_loadCh(basepath, basename, session, edCh, win, bit2uv);
end

% the 'eeg' branch reads the whole trace; window it to match 'ripp', which
% evt_loadCh already windows on the way in
if strcmp(chMode, 'eeg')
    sig = winCrop(sig, win, fs);
end

%% ========================================================================
%  EMG  (windowed to the detection signal; dropped on length mismatch)
%  ========================================================================
emg = [];
if isfile(sigFile)
    S = load(sigFile, 'emg');
    if isfield(S, 'emg'), emg = S.emg(:); end
end

if ~isempty(emg)
    emg = winCrop(emg, win, fs);
    if numel(emg) ~= numel(sig)
        warning('ed_sigLoad:emgFit', ...
            ['EMG (%d) and signal (%d) lengths differ; ', ...
            'skipping the EMG metric.'], ...
            numel(emg), numel(sig));
        emg = [];
    end
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function y = winCrop(x, win, fs)
% Crop a full-recording trace to the analysis window at rate fs.
s1 = max(1, round(win(1) * fs) + 1);
if isinf(win(2))
    s2 = numel(x);
else
    s2 = min(numel(x), round(win(2) * fs));
end
y = x(s1 : s2);
end
