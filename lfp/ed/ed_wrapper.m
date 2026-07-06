function ed = ed_wrapper(varargin)
% ED_WRAPPER Detect, characterise, and curate epileptiform discharges (EDs).
%
%   ed = ED_WRAPPER(varargin)
%
%   SUMMARY:
%       Orchestrates the ED pipeline for one session (mirrors ripp_wrapper,
%       struct-based, no class). The body reads top-to-bottom:
%       1. Setup: load session + sleep states; prepare spike times (evt_spkPrep)
%          and bout times (evt_boutTimes).
%       2. Detect-or-load: if <basename>.ed.mat exists and flgForce is false,
%          load it (the cheap re-curate path); otherwise:
%          a. Signal: load the detection signal + EMG (ed_sigLoad).
%          b. Detect (ed_detect); per-event features (ed_params); QA on EMG
%             (+ optional amp / dur; no state criterion) via evt_qa as a filter
%             (evt_subset drops the failures).
%          c. Characterise: matched controls (evt_ctrlTimes), vigilance state
%             (evt_states), LFP maps (evt_maps).
%          d. Spiking: SU/MU modulation + PETH (evt_spks).
%          e. Convert to absolute time, populate .info, save, plot.
%       3. Optionally launch the curation GUI (guiPath_curate, preset 'EDs').
%
%   INPUTS (Parameter/Value):
%       'basepath'   - (Char) Session directory. {pwd}
%       'win'        - (Vec)  Analysis window [start end] (s). {[0 Inf]}
%       'sigSource'  - (Char) 'lfp' (raw channel, default) | 'eeg' (sSig.eeg).
%       'edCh'       - (Num)  Channel for 'lfp'. {channelTags.Ripple}
%       'bit2uv'     - (Num)  Conversion for the 'lfp' source. {auto}
%       'thr'        - (Num)  Z-score detection threshold. {7}
%       'thrDir'     - (Char) 'positive' | 'negative' | 'both'. {'both'}
%       'baseWin'    - (Num)  Moving baseline window [s]. {5}
%       'interDur'   - (Num)  Refractory / burst-merge window [s]. {0.025}
%       'ampWin'     - (Num)  Peak-to-peak gate half-window [s]. {0.015}
%       'lowThr'     - (Num)  Twin-peak merge trough (signal units). {0.2}
%       'minAmp'     - (Num)  Absolute amplitude floor at detection. {[]}
%       'marg'       - (Num)  Feature clip half-window [s]. {0.05}
%       'thrZ'       - (Num)  EMG z pass threshold for QA. {3}
%       'durLim'     - (Vec)  [min max] half-amp width for QA (ms). {[]}
%       'minAmpQA'   - (Num)  Amplitude floor for QA. {[]}
%       'mapDur'     - (Vec)  PETH / map window [pre post] (s). {[-0.1 0.1]}
%       'flgPlot'    - (Log)  Generate the summary figure? {true}
%       'flgSave'    - (Log)  Save output .mat files? {false}
%       'flgCurate'  - (Log)  Launch the curation GUI? {false}
%       'flgForce'   - (Log)  Re-detect even if .ed.mat exists? {false}
%       'verbose'    - (Log)  Print progress? {true}
%
%   OUTPUT:
%       ed           - (Struct) Events + per-event metrics (canonical schema):
%           .times .peakTime .pos .state .accepted .ctrlTimes .info, discharge
%           params (.amp .ampZ .dur .width10), QA metric (.emgZ), and population
%           spike metrics (.spks). When the GUI is launched, curated acceptance
%           is saved by the GUI; this return is the pre-curation struct.
%
%   FILES SAVED (when flgSave = true):
%       basename.ed.mat         - events + per-event population metrics
%       basename.edStates.mat   - per-bout rate/density table (evt_states)
%       basename.edMaps.mat     - per-event signal maps
%       basename.edSpks.mat     - per-unit spike stats + PETH (light)
%       basename.edSpkMaps.mat  - 3D spike raster [unit x event x bin] (heavy)
%
%   DEPENDENCIES:
%       ed_sigLoad, ed_detect, ed_params, guiPath_curate, basepaths2vars; and
%       the shared event layer (lfp/events): evt_files, evt_pickCh, evt_spkPrep,
%       evt_boutTimes, evt_emgScore, evt_qa, evt_subset, evt_ctrlTimes,
%       evt_states, evt_maps, evt_spks, evt_saveSpks, evt_plotSpks.
%
%   HISTORY:
%       Created: 260622
%       Updated: 260706 (shared spine via evt_* helpers; expose mapDur; signal
%                load moved into the detect branch).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'sigSource', 'lfp', @(x) any(strcmpi(x, {'eeg', 'lfp'})));
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'thr', 7, @isnumeric);
addParameter(p, 'thrDir', 'both', @ischar);
addParameter(p, 'baseWin', 5, @isnumeric);
addParameter(p, 'interDur', 0.025, @isnumeric);
addParameter(p, 'ampWin', 0.015, @isnumeric);
addParameter(p, 'lowThr', 0.2, @isnumeric);
addParameter(p, 'minAmp', [], @isnumeric);
addParameter(p, 'marg', 0.05, @isnumeric);
addParameter(p, 'thrZ', 3, @isnumeric);
addParameter(p, 'durLim', [], @isnumeric);
addParameter(p, 'minAmpQA', [], @isnumeric);
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

basepath  = p.Results.basepath;
win       = p.Results.win;
sigSource = lower(p.Results.sigSource);
edCh      = p.Results.edCh;
bit2uv    = p.Results.bit2uv;
thr       = p.Results.thr;
thrDir    = p.Results.thrDir;
baseWin   = p.Results.baseWin;
interDur  = p.Results.interDur;
ampWin    = p.Results.ampWin;
lowThr    = p.Results.lowThr;
minAmp    = p.Results.minAmp;
marg      = p.Results.marg;
thrZ      = p.Results.thrZ;
durLim    = p.Results.durLim;
minAmpQA  = p.Results.minAmpQA;
mapDur    = p.Results.mapDur;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgCurate = p.Results.flgCurate;
flgForce  = p.Results.flgForce;
verbose   = p.Results.verbose;

%% ========================================================================
%  SETUP
%  ========================================================================
basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
files = evt_files(basepath, basename, 'ed');

if verbose, fprintf('[ED]: Session %s\n', basename); end

% Session + data. Session may be absent for minimal (eeg-only) layouts, so it
% is guarded here and everywhere it is read.
v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'sleep_states', 'spikes', 'spktimes', 'units'});
session = [];
if isfield(v, 'session'), session = v.session; end
if isinf(win(2)), sigDur = Inf; else, sigDur = win(2) - win(1); end

% Spike times (window-relative SU + pooled MUA). fsSpk is the wideband rate
% (spktimes are wideband samples); absent when the session lacks it, in which
% case the spike analyses are skipped.
fsSpk = [];
if ~isempty(session) && isfield(session, 'extracellular') ...
        && isfield(session.extracellular, 'sr')
    fsSpk = session.extracellular.sr;
end
spkTimes = {}; muTimes = {[]}; uType = [];
if ~isempty(fsSpk)
    [spkTimes, muTimes, uType] = evt_spkPrep(v, win, sigDur, fsSpk);
end
hasSpikes = ~isempty(spkTimes);

% Bout times (window-relative). ED ignores the NREM baseline (3rd output).
[boutTimes, vldTimes, ~] = evt_boutTimes(v, win, sigDur);

%% ========================================================================
%  DETECT OR LOAD
%  ========================================================================
if isfile(files.evt) && ~flgForce
    if verbose, fprintf('[ED]: Loading existing %s.ed.mat\n', basename); end
    S = load(files.evt, 'ed');
    ed = S.ed;

else
    % ---- Signal ---------------------------------------------------------
    % Detect on the ripple-tagged channel so the ED signal is the same LFP as
    % the ripples pipeline (only for the 'lfp' source; 'eeg' ignores edCh).
    if strcmp(sigSource, 'lfp') && isempty(edCh)
        edCh = evt_pickCh(session);
    end
    if verbose, fprintf('[ED]: Loading signal...\n'); end
    [sig, emg, ~, fs] = ed_sigLoad(basepath, 'sigSource', sigSource, ...
        'win', win, 'session', session, 'basename', basename, ...
        'edCh', edCh, 'bit2uv', bit2uv);

    % ---- Detect + QA ----------------------------------------------------
    if verbose, fprintf('[ED]: Detecting...\n'); end
    ed = ed_detect(sig, fs, 'thr', thr, 'thrDir', thrDir, ...
        'baseWin', baseWin, 'interDur', interDur, 'ampWin', ampWin, ...
        'lowThr', lowThr, 'minAmp', minAmp);
    ed = ed_params(sig, ed, 'marg', marg);

    % QA metric: EMG (mean event EMG vs the valid-state baseline; whole
    % recording when states are absent). Amplitude / duration gates optional.
    edWin  = 0.05;   % EMG measurement half-window around each peak [s]
    edWins = [ed.peakTime(:) - edWin, ed.peakTime(:) + edWin];
    ed.emgZ = evt_emgScore(emg, edWins, fs, 'baselineTimes', vldTimes);

    metrics = ed.emgZ;
    ranges  = {[-Inf, thrZ]};
    names   = {'emg'};
    if ~isempty(durLim)
        metrics(:, end+1) = ed.dur(:);
        ranges{end+1} = durLim;              % [min max]
        names{end+1}  = 'dur';
    end
    if ~isempty(minAmpQA)
        metrics(:, end+1) = ed.amp(:);
        ranges{end+1} = [minAmpQA, Inf];
        names{end+1}  = 'amp';
    end

    % QA filter (no state criterion for EDs). Drop the failures, remove any
    % transient .idxQA, and seed acceptance all-true on the survivors.
    keep = evt_qa(ed.peakTime, 'metrics', metrics, ...
        'ranges', ranges, 'names', names);
    ed = evt_subset(ed, keep);
    if isfield(ed, 'idxQA'), ed = rmfield(ed, 'idxQA'); end
    ed.accepted = true(numel(ed.pos), 1);

    % ---- Characterise ---------------------------------------------------
    ed.ctrlTimes = evt_ctrlTimes(ed.times, ...
        'vldTimes', vldTimes, 'flgPlot', false);

    nEvt = numel(ed.pos);
    if ~isempty(boutTimes) && nEvt > 0
        [ed.state, ~] = evt_states(ed.times, ed.peakTime, boutTimes, ...
            'basepath', basepath, 'flgSave', flgSave, 'flgPlot', false, ...
            'name', 'ed', 'lbl', 'ED');
    else
        ed.state = categorical(nan(nEvt, 1));
    end

    if verbose, fprintf('[ED]: Generating signal maps...\n'); end
    edMaps = evt_maps(struct('lfp', sig(:)), ed.peakTime, fs, ...
        'mapDur', mapDur, 'flgSave', false);

    % ---- Spiking (SU/MU modulation + PETH) ------------------------------
    hasSpks = hasSpikes && nEvt > 0;
    if hasSpks
        if verbose, fprintf('[ED]: Analysing spikes...\n'); end
        edSpks = evt_spks(spkTimes, muTimes, ed.times, ed.ctrlTimes, ...
            ed.peakTime, 'unitType', uType, 'mapDur', mapDur, 'winFxd', 0.050);
        ed.spks = edSpks.events;              % per-event metrics live on ed
        edSpks = rmfield(edSpks, 'events');
    end

    % ---- Finalise: absolute time + provenance ---------------------------
    ed.times     = ed.times + win(1);
    ed.peakTime  = ed.peakTime + win(1);
    ed.pos       = ed.pos + round(win(1) * fs);
    ed.ctrlTimes = ed.ctrlTimes + win(1);

    ed.info.basename  = basename;
    ed.info.sigSource = sigSource;
    ed.info.edCh      = edCh;        % detection channel; the GUI matches it
    ed.info.win       = win;
    ed.info.runtime   = datetime('now');

    % ---- Save -----------------------------------------------------------
    if flgSave
        if verbose, fprintf('[ED]: Saving output files...\n'); end
        save(files.evt, 'ed', '-v7.3');
        save(files.maps, 'edMaps', '-v7.3');
        if hasSpks
            evt_saveSpks(edSpks, basepath, basename, 'ed');
        end
    end

    % ---- Plot -----------------------------------------------------------
    if flgPlot && hasSpks
        evt_plotSpks(edSpks, 'basepath', basepath, 'flgSaveFig', true, ...
            'name', 'ed', 'lbl', 'ED');
    end
end

%% ========================================================================
%  CURATE
%  ========================================================================
if flgCurate
    if isfile(files.evt)
        if verbose
            fprintf('[ED]: Launching curation GUI (%d events)\n', ...
                numel(ed.pos));
        end
        guiPath_curate(basepath, 'preset', 'EDs', 'basename', basename);
    elseif verbose
        fprintf('[ED]: No %s.ed.mat on disk; set flgSave=true to curate.\n', ...
            basename);
    end
end

if verbose, fprintf('[ED]: Done (%s).\n', basename); end

end     % EOF
