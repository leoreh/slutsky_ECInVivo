function ed = ed_wrapper(varargin)
% ED_WRAPPER Detect, characterise, and curate epileptiform discharges (EDs).
%
%   ed = ED_WRAPPER(varargin)
%
%   SUMMARY:
%       Orchestrates the ED pipeline for one recording session (mirrors the
%       ripples pipeline, struct-based, no class):
%       1. Load session + sleep_states; load signals once (ed_sigLoad).
%       2. Detect (ed_detect), characterise (ed_params), score EMG
%          (ed_reject_emg), label state (ed_states).
%       3. Seed acceptance from quality masks (idxQA); never delete events.
%       4. Convert times to absolute, optionally save <basename>.ed.mat.
%       5. Optionally launch the curation GUI (ed_gui).
%       If <basename>.ed.mat already exists and flgForce is false, detection
%       is skipped and the stored result is loaded straight into the GUI
%       (the cheap re-curate path).
%
%   INPUTS (Parameter/Value):
%       'basepath'   - (Char) Session directory. {pwd}
%       'win'        - (Vec)  Analysis window [start end] (s). {[0 Inf]}
%       'sigSource'  - (Char) 'eeg' (sSig.eeg, default) | 'lfp' (raw channel).
%       'edCh'       - (Num)  Channel for the 'lfp' source. {1}
%       'bit2uv'     - (Num)  Conversion for the 'lfp' source. {auto}
%       'thr'        - (Num)  Z-score detection threshold. {7}
%       'thrDir'     - (Char) 'positive' | 'negative' | 'both'. {'both'}
%       'baseWin'    - (Num)  Moving baseline window [s]. {5}
%       'interDur'   - (Num)  Refractory / burst-merge window [s]. {0.025}
%       'ampWin'     - (Num)  Peak-to-peak gate half-window [s]. {0.015}
%       'lowThr'     - (Num)  Twin-peak merge trough (signal units). {0.2}
%       'minAmp'     - (Num)  Absolute amplitude floor at detection. {[]}
%       'marg'       - (Num)  Feature clip half-window [s]. {0.05}
%       'emgMethod'  - (Char) 'zscore' (default) | 'emg_rms'.
%       'thrZ'       - (Num)  EMG z-score pass threshold. {3}
%       'thrRms'     - (Num)  emg_rms pass threshold. {75th pct}
%       'durLim'     - (Vec)  [min max] half-amp width for QA (ms). {[]}
%       'minAmpQA'   - (Num)  Amplitude floor for QA. {[]}
%       'binsize'    - (Num)  Bin width for ed_rate [s]. {60}
%       'flgRate'    - (Log)  Compute (and, with flgPlot, plot) ed_rate? {false}
%       'flgPlot'    - (Log)  Launch the curation GUI? {true}
%       'flgSave'    - (Log)  Save <basename>.ed.mat (+ edStates)? {false}
%       'flgForce'   - (Log)  Re-detect even if <basename>.ed.mat exists? {false}
%       'verbose'    - (Log)  Print progress? {true}
%
%   OUTPUT:
%       ed          - (Struct) Detection result (see ed_params / schema). When
%                     the GUI is launched, curated acceptance is saved to disk
%                     by the GUI; this return is the pre-curation struct.
%
%   DEPENDENCIES:
%       basepaths2vars, ed_sigLoad, ed_detect, ed_params, ed_reject_emg,
%       ed_states, ed_rate, ed_gui.
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'sigSource', 'eeg', @(x) any(strcmpi(x, {'eeg', 'lfp'})));
addParameter(p, 'edCh', 1, @isnumeric);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'thr', 7, @isnumeric);
addParameter(p, 'thrDir', 'both', @ischar);
addParameter(p, 'baseWin', 5, @isnumeric);
addParameter(p, 'interDur', 0.025, @isnumeric);
addParameter(p, 'ampWin', 0.015, @isnumeric);
addParameter(p, 'lowThr', 0.2, @isnumeric);
addParameter(p, 'minAmp', [], @isnumeric);
addParameter(p, 'marg', 0.05, @isnumeric);
addParameter(p, 'emgMethod', 'zscore', @ischar);
addParameter(p, 'thrZ', 3, @isnumeric);
addParameter(p, 'thrRms', [], @isnumeric);
addParameter(p, 'durLim', [], @isnumeric);
addParameter(p, 'minAmpQA', [], @isnumeric);
addParameter(p, 'binsize', 60, @isnumeric);
addParameter(p, 'flgRate', false, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);

parse(p, varargin{:});
basepath  = p.Results.basepath;
win       = p.Results.win;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgForce  = p.Results.flgForce;
flgRate   = p.Results.flgRate;
verbose   = p.Results.verbose;

%% ========================================================================
%  SETUP
%  ========================================================================
cd(basepath);
basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
edFile = fullfile(basepath, [basename, '.ed.mat']);

if verbose, fprintf('[ED]: Session %s\n', basename); end

% Session metadata + sleep states (absent for some layouts; handled below)
v = basepaths2vars('basepaths', {basepath}, 'vars', {'session', 'sleep_states'});
session = [];
if isfield(v, 'session'), session = v.session; end

% Signals (loaded once; sSig/specAdapter are full-session for the GUI)
[sig, emg, emgRms, fs, specAdapter, sSig] = ed_sigLoad(basepath, ...
    'sigSource', p.Results.sigSource, 'win', win, 'session', session, ...
    'basename', basename, 'edCh', p.Results.edCh, 'bit2uv', p.Results.bit2uv);

% Duration for bout clipping
if isinf(win(2)), sigDur = Inf; else, sigDur = win(2) - win(1); end

% Bout times (relative to win), for state assignment
boutTimes = [];
if isfield(v, 'ss') && isfield(v.ss, 'bouts') && isfield(v.ss.bouts, 'times')
    boutTimes = v.ss.bouts.times;
    boutTimes = cellfun(@(x) x - win(1), boutTimes, 'uni', false);
    boutTimes = cellfun(@(x) x(x(:,2) > 0 & x(:,1) < sigDur, :), boutTimes, 'uni', false);
else
    warning('ed_wrapper:noStates', 'sleep_states not found; events left unlabelled.');
end

%% ========================================================================
%  DETECT OR LOAD
%  ========================================================================

if isfile(edFile) && ~flgForce
    if verbose, fprintf('[ED]: Loading existing %s (skip detection)\n', [basename, '.ed.mat']); end
    S = load(edFile, 'ed');
    ed = S.ed;

else
    if verbose, fprintf('[ED]: Detecting...\n'); end
    ed = ed_detect(sig, fs, ...
        'thr', p.Results.thr, 'thrDir', p.Results.thrDir, ...
        'baseWin', p.Results.baseWin, 'interDur', p.Results.interDur, ...
        'ampWin', p.Results.ampWin, 'lowThr', p.Results.lowThr, ...
        'minAmp', p.Results.minAmp);

    % Per-event features
    ed = ed_params(sig, ed, 'marg', p.Results.marg);

    % EMG scoring
    if strcmpi(p.Results.emgMethod, 'zscore')
        ed = ed_reject_emg(ed, emg, fs, 'method', 'zscore', 'thrZ', p.Results.thrZ);
    else
        ed = ed_reject_emg(ed, [], fs, 'method', 'emg_rms', ...
            'emgRms', emgRms, 'thrRms', p.Results.thrRms);
    end

    % State assignment (relative frame, before absolute conversion)
    nEvt = numel(ed.pos);
    if ~isempty(boutTimes) && nEvt > 0
        [ed.state, ~] = ed_states(ed, boutTimes, 'basepath', basepath, 'flgSave', flgSave);
    else
        ed.state = categorical(nan(nEvt, 1));
    end

    % Quality masks -> seed acceptance (never delete)
    if isempty(p.Results.durLim)
        ed.idxQA.dur = true(nEvt, 1);
    else
        ed.idxQA.dur = ed.dur >= p.Results.durLim(1) & ed.dur <= p.Results.durLim(2);
    end
    if isempty(p.Results.minAmpQA)
        ed.idxQA.amp = true(nEvt, 1);
    else
        ed.idxQA.amp = ed.amp >= p.Results.minAmpQA;
    end
    ed.idxQA.auto = ed.idxQA.emg & ed.idxQA.amp & ed.idxQA.dur;
    ed.accepted = ed.idxQA.auto(:);

    % Absolute times (events live in the full-session frame the GUI shows)
    ed.times    = ed.times + win(1);
    ed.peakTime = ed.peakTime + win(1);
    ed.pos      = ed.pos + round(win(1) * fs);

    % Finalise info
    ed.info.basename  = basename;
    ed.info.sigSource = p.Results.sigSource;
    ed.info.win       = win;
    ed.info.runtime   = datetime('now');

    % Optional rate
    if flgRate
        ed.rate = ed_rate(ed, 'binsize', p.Results.binsize, ...
            'flgPlot', flgPlot, 'basepath', basepath);
    end

    % Save
    if flgSave
        if verbose, fprintf('[ED]: Saving %s\n', [basename, '.ed.mat']); end
        save(edFile, 'ed', '-v7.3');
    end
end

%% ========================================================================
%  CURATION GUI
%  ========================================================================

if flgPlot
    if verbose, fprintf('[ED]: Launching curation GUI (%d events)\n', numel(ed.pos)); end
    ed_gui(basepath, ed, sSig, 'specAdapter', specAdapter, 'basename', basename);
end

if verbose, fprintf('[ED]: Done (%s).\n', basename); end

end     % EOF
