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
%          (ed_reject_emg), label state (evt_states).
%       3. Filter by automatic QA (EMG + optional amp/dur); seed .accepted
%          all-true on the survivors (no .idxQA kept).
%       4. Parity analyses (shared evt_* layer): matched control intervals,
%          LFP maps, and MUA/SU spike modulation + PETH around discharges.
%       5. Convert times to absolute; optionally save .ed / .edMaps /
%          .edSpks (per-unit stats + 3D PETH, consolidated).
%       6. Optionally launch the curation GUI (gui_curate, preset 'EDs').
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
%       'binsize'    - (Num)  Bin width for evt_rate [s]. {60}
%       'flgRate'    - (Log)  Compute (and, with flgPlot, plot) evt_rate? {false}
%       'flgPlot'    - (Log)  Generate summary figures + viewers? {true}
%       'flgSave'    - (Log)  Save <basename>.ed.mat (+ artefacts)? {false}
%       'flgCurate'  - (Log)  Launch the curation GUI (gui_curate)? {false}
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
%       evt_spkPrep, evt_states, evt_ctrlTimes, evt_maps, evt_spkAnalysis,
%       evt_plotSpks, evt_viewSpks, evt_rate, gui_curate.
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
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);

parse(p, varargin{:});
basepath  = p.Results.basepath;
win       = p.Results.win;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgCurate = p.Results.flgCurate;
flgForce  = p.Results.flgForce;
flgRate   = p.Results.flgRate;
verbose   = p.Results.verbose;

%% ========================================================================
%  SETUP
%  ========================================================================
basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end
edFile = fullfile(basepath, [basename, '.ed.mat']);

if verbose, fprintf('[ED]: Session %s\n', basename); end

% Session metadata + sleep states (absent for some layouts; handled below)
v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'sleep_states', 'spikes', 'spktimes', 'units'});
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

% Valid-state intervals for control matching (QWAKE+LSLEEP+NREM, mirroring
% ripples). Empty when states are absent -> controls span the whole recording.
vldTimes = [];
if ~isempty(boutTimes) && numel(boutTimes) >= 4
    vldTimes = vertcat(boutTimes{2}, boutTimes{3}, boutTimes{4});
end

% Spikes (optional) for the parity analyses, prepared by the shared layer
% (identical treatment to the ripple pipeline): window-relative single-unit
% times, one pooled MUA vector, and unit types. Absent on sessions without
% sorted spikes, in which case the spike analyses are skipped.
fsSpk = fs;
if ~isempty(session) && isfield(session, 'extracellular') && isfield(session.extracellular, 'sr')
    fsSpk = session.extracellular.sr;
end
spikesIn = []; if isfield(v, 'spikes'), spikesIn = v.spikes; end
spktimesIn = []; if isfield(v, 'spktimes'), spktimesIn = v.spktimes; end
unitsIn = []; if isfield(v, 'units'), unitsIn = v.units; end
[spkTimes, muTimes, uType] = evt_spkPrep(spikesIn, spktimesIn, unitsIn, ...
    win, sigDur, fsSpk);
hasSpikes = ~isempty(spkTimes);

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
        [ed.state, ~] = evt_states(ed.times, ed.peakTime, boutTimes, ...
            'basepath', basepath, 'flgSave', flgSave, 'flgPlot', false, 'name', 'ed');
    else
        ed.state = categorical(nan(nEvt, 1));
    end

    % Automatic QA as a FILTER: build the pass mask (EMG from ed_reject_emg,
    % plus optional amplitude / duration gates), drop the failures, and seed
    % acceptance all-true on the survivors. The per-criterion breakdown is
    % applied then discarded (no .idxQA stored) - mirrors the ripple pipeline
    % so only accepted candidates plus the curation mask persist on disk.
    qaPass = ed.idxQA.emg(:);
    if ~isempty(p.Results.durLim)
        qaPass = qaPass & ed.dur(:) >= p.Results.durLim(1) ...
                        & ed.dur(:) <= p.Results.durLim(2);
    end
    if ~isempty(p.Results.minAmpQA)
        qaPass = qaPass & ed.amp(:) >= p.Results.minAmpQA;
    end
    ed = ed_filterEvents(ed, qaPass);
    ed.accepted = true(numel(ed.pos), 1);

    % ==== Parity analyses (relative frame, mirroring the ripple pipeline) ====
    mapDur = [-0.1 0.1];

    % Matched control intervals from valid vigilance states
    ed.ctrlTimes = evt_ctrlTimes(ed.times, 'vldTimes', vldTimes, 'flgPlot', false);

    % LFP maps around each discharge peak (detection signal)
    edMaps = evt_maps(struct('lfp', sig(:)), ed.peakTime, fs, ...
        'mapDur', mapDur, 'flgSave', false);

    % MUA / single-unit spike modulation + PETH (only when spikes exist),
    % consolidated into one struct by the shared layer.
    edSpks = struct(); flgEdSpks = false;
    if hasSpikes && numel(ed.pos) > 0
        flgEdSpks = true;
        edSpks = evt_spkAnalysis(spkTimes, muTimes, ed.times, ed.ctrlTimes, ...
            ed.peakTime, 'unitType', uType, 'mapDur', mapDur, 'winFxd', 0.050);

        % Move per-discharge population metrics onto the ed struct
        ed.spks = edSpks.events;
        edSpks = rmfield(edSpks, 'events');
    end

    % Absolute times (events live in the full-session frame the GUI shows)
    ed.times     = ed.times + win(1);
    ed.peakTime  = ed.peakTime + win(1);
    ed.pos       = ed.pos + round(win(1) * fs);
    ed.ctrlTimes = ed.ctrlTimes + win(1);

    % Finalise info
    ed.info.basename  = basename;
    ed.info.sigSource = p.Results.sigSource;
    ed.info.win       = win;
    ed.info.runtime   = datetime('now');

    % Optional rate
    if flgRate
        ed.rate = evt_rate(ed, 'binsize', p.Results.binsize, ...
            'flgPlot', flgPlot, 'basepath', basepath);
    end

    % Save (event struct + parity artefacts)
    if flgSave
        if verbose, fprintf('[ED]: Saving %s\n', [basename, '.ed.mat']); end
        save(edFile, 'ed', '-v7.3');
        save(fullfile(basepath, [basename, '.edMaps.mat']), 'edMaps', '-v7.3');
        if flgEdSpks
            save(fullfile(basepath, [basename, '.edSpks.mat']), 'edSpks', '-v7.3');
        end
    end

    % Spike-modulation summary figure + interactive viewers (parity with
    % ripples, from the shared layer)
    if flgPlot && flgEdSpks
        evt_plotSpks(edSpks, 'basepath', basepath, ...
            'flgSaveFig', true, 'name', 'ed', 'lbl', 'ED');

        edState = [];
        if isfield(ed, 'state'), edState = ed.state; end
        evt_viewSpks(edMaps, edSpks, uType, edState, 'mapYVar', 'lfp');
    end
end

%% ========================================================================
%  CURATION GUI
%  ========================================================================

if flgCurate
    if isfile(edFile)
        if verbose, fprintf('[ED]: Launching curation GUI (%d events)\n', numel(ed.pos)); end
        gui_curate(basepath, 'preset', 'EDs', 'basename', basename);
    elseif verbose
        fprintf('[ED]: No %s on disk; set flgSave=true to curate.\n', [basename, '.ed.mat']);
    end
end

if verbose, fprintf('[ED]: Done (%s).\n', basename); end

end     % EOF


%% ========================================================================
%  HELPER: ED_FILTEREVENTS
%  ========================================================================
function ed = ed_filterEvents(ed, keep)
% Subset every per-event field of the ed struct by the logical KEEP mask and
% drop the transient .idxQA breakdown. Per-event fields are those whose first
% dimension equals the event count; scalar / struct fields (.info) are left.
keep = logical(keep(:));
nEvt = numel(keep);
fn = fieldnames(ed);
for iF = 1:numel(fn)
    f = ed.(fn{iF});
    if (isnumeric(f) || islogical(f) || iscategorical(f)) && size(f, 1) == nEvt
        ed.(fn{iF}) = f(keep, :);
    end
end
if isfield(ed, 'idxQA'), ed = rmfield(ed, 'idxQA'); end
end     % EOF
