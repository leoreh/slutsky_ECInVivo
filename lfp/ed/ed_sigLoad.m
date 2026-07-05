function [sig, emg, emgRms, fs, specAdapter, sSig] = ed_sigLoad(basepath, varargin)
% ED_SIGLOAD Load the signals for ED detection / curation, once.
%
%   [sig, emg, emgRms, fs, specAdapter, sSig] = ED_SIGLOAD(basepath, varargin)
%
%   SUMMARY:
%       Loads <basename>.sleep_sig.mat a single time and returns the pieces
%       the ED pipeline and GUI need. The detection signal (SIG) and EMG are
%       windowed to 'win'; SSIG and SPECADAPTER are returned whole-session so
%       the GUI can show full-session context against absolute event times.
%       The default detection signal is sSig.eeg; an 'lfp' escape hatch loads
%       a raw channel via binary_load.
%
%   INPUTS:
%       basepath    - (Char) Session directory holding <basename>.sleep_sig.mat.
%       varargin    - Parameter/Value pairs:
%           'sigSource' - (Char) 'eeg' (default) | 'lfp'.
%           'win'       - (Vec)  Window [start end] (s) for SIG/EMG. {[0 Inf]}
%           'session'   - (Struct) Session metadata (loaded if empty & needed).
%           'edCh'      - (Num)  Channel for the 'lfp' source. {1}
%           'bit2uv'    - (Num)  Conversion for the 'lfp' source. {auto}
%
%   OUTPUTS:
%       sig         - (Vec)    Windowed detection signal.
%       emg         - (Vec)    Windowed EMG (matches SIG length).
%       emgRms      - (Vec)    Full-session 1-Hz log-RMS EMG (or []).
%       fs          - (Num)    Sampling frequency of SIG/EMG [Hz].
%       specAdapter - (Struct) Full-session spectrogram for plot_spec
%                              (.s/.freq/.tstamps).
%       sSig        - (Struct) The full sleep_sig struct (for the GUI).
%
%   DEPENDENCIES:
%       binary_load, basepaths2vars (only for the 'lfp' source).
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'sigSource', 'eeg', @(x) any(strcmpi(x, {'eeg', 'lfp'})));
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'edCh', 1, @isnumeric);
addParameter(p, 'bit2uv', [], @isnumeric);

parse(p, basepath, varargin{:});
basepath  = p.Results.basepath;
sigSource = lower(p.Results.sigSource);
win       = p.Results.win;
session   = p.Results.session;

% basename defaults to the folder name (standard one-session-per-folder
% layout) but can be overridden when files are named differently from the
% containing folder (e.g. several recordings stored flat in one directory).
basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
if ~isfile(sigfile)
    error('ed_sigLoad:noFile', 'sleep_sig.mat not found: %s', sigfile);
end

%% ========================================================================
%  LOAD SSIG (ONCE)
%  ========================================================================

% sleep_sig.mat is saved with -struct, so its fields are top-level variables
sSig   = load(sigfile);
fs     = sSig.fs;
emgAll = sSig.emg(:);

if isfield(sSig, 'emg_rms')
    emgRms = sSig.emg_rms(:);
else
    emgRms = [];
end

% Spectrogram adapter for plot_spec (references, not copies)
specAdapter = struct('s', sSig.spec, 'freq', sSig.spec_freq(:), ...
    'tstamps', sSig.spec_tstamps(:));

%% ========================================================================
%  DETECTION SIGNAL
%  ========================================================================

switch sigSource
    case 'eeg'
        sigAll = sSig.eeg(:);

    case 'lfp'
        % Escape hatch: detect on a raw LFP channel instead of sSig.eeg.
        % binary_load applies bit2uv once, so sigAll is in uV; its absolute
        % scale may differ from sSig.eeg, but detection is z-scored so this
        % affects only reported amplitude units, not which events are found.
        if isempty(session)
            v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'});
            session = v.session;
        end
        fsLfp = session.extracellular.srLfp;
        if abs(fsLfp - fs) > 1
            warning('ed_sigLoad:fsMismatch', ...
                'srLfp (%g) differs from sSig.fs (%g); EMG/spec may misalign.', fsLfp, fs);
        end
        fs = fsLfp;
        nChans = session.extracellular.nChannels;
        bit2uv = p.Results.bit2uv;
        if isempty(bit2uv)
            if round(session.extracellular.sr) == 24414
                bit2uv = 1;        % TDT
            else
                bit2uv = 0.195;    % Intan
            end
        end
        fname = fullfile(basepath, [basename, '.lfp']);
        sigAll = double(binary_load(fname, 'duration', Inf, 'fs', fs, ...
            'nCh', nChans, 'start', 0, 'ch', p.Results.edCh, ...
            'downsample', 1, 'bit2uv', bit2uv));
        if size(sigAll, 2) > 1
            sigAll = mean(sigAll, 2);
        end
        sigAll = sigAll(:);
end

%% ========================================================================
%  WINDOW SIG / EMG
%  ========================================================================

s1 = max(1, round(win(1) * fs) + 1);
if isinf(win(2))
    s2sig = numel(sigAll);
    s2emg = numel(emgAll);
else
    s2sig = min(numel(sigAll), round(win(2) * fs));
    s2emg = min(numel(emgAll), round(win(2) * fs));
end
sig = sigAll(s1:s2sig);
emg = emgAll(s1:s2emg);

end     % EOF
