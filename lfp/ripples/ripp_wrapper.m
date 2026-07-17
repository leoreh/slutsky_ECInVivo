function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Detect, characterise, and curate sharp-wave ripples (SWR).
%
%   ripp = RIPP_WRAPPER(varargin)
%
%   SUMMARY:
%       Orchestrates the ripple pipeline for one session (mirrors ed_wrapper,
%       struct-based, no class). The body reads top-to-bottom:
%       1. Setup: load session + sleep states; prepare spike times (evt_spkPrep)
%          and bout times (evt_boutTimes).
%       2. Detect-or-load: if <basename>.ripp.mat exists and flgForce is false,
%          load it (the cheap re-curate path); otherwise:
%          a. Signal: load the ripple channel + EMG (ripp_sigLoad); filter +
%             envelope + z-score (ripp_sigPrep).
%          b. Detect (ripp_times); QA on state / EMG / spike-gain (evt_qa) as a
%             filter (evt_subset drops the failures).
%          c. Characterise: matched controls (evt_ctrlTimes), vigilance state
%             (evt_states), per-event params (ripp_params), LFP maps (evt_maps).
%          d. Spiking: SU/MU modulation + PETH (evt_spks); phase (spklfp_phase).
%          e. Convert to absolute time, populate .info, save, plot.
%       3. Optionally export a NeuroScope event file (evt2ns). This sits
%          outside the detect-or-load branch, so re-running on a curated
%          session refreshes the export against the stored .accepted mask.
%       4. Optionally launch the curation GUI (guiPath, 'Ripples').
%
%   INPUTS (Parameter/Value):
%       'basepath'   - (Char) Session directory. {pwd}
%       'rippCh'     - (Num)  Zero-indexed channel. {channelTags.Ripple}
%       'thr'        - (Vec)  Thresholds [start peak cont max min_cont].
%       'limDur'     - (Vec)  Duration limits [min max inter min_cont] (ms).
%       'thrEmg'     - (Num)  EMG z pass threshold for QA. {2}
%       'win'        - (Vec)  Analysis window [start end] (s). {[0 Inf]}
%       'passband'   - (Vec)  Filter band [min max] (Hz). {[80 250]}
%       'detectMet'  - (Num)  Detection method id (see ripp_sigPrep). {3}
%       'zMet'       - (Char) Z-scoring method ('adaptive', 'nrem'). {'nrem'}
%       'mapDur'     - (Vec)  PETH / map window [pre post] (s). {[-0.1 0.1]}
%       'bit2uv'     - (Num)  Conversion factor. {auto: TDT vs Intan}
%       'flgPlot'    - (Log)  Generate the summary figure? {true}
%       'flgSave'    - (Log)  Save output .mat files? {false}
%       'flgNS'      - (Log)  Write the NeuroScope event file? {false}
%       'flgCurate'  - (Log)  Launch the curation GUI? {false}
%       'flgForce'   - (Log)  Re-detect even if .ripp.mat exists? {false}
%       'verbose'    - (Log)  Print progress? {true}
%
%   OUTPUT:
%       ripp         - (Struct) Events + per-event metrics (canonical schema):
%           .times .peakTime .state .accepted .ctrlTimes .info, ripple params
%           (.amp .freq .freqEvent .energy .dur .skew), QA metrics (.emg
%           .spkGain), and population spike metrics (.spks).
%
%   FILES SAVED (when flgSave = true):
%       basename.ripp.mat        - events + per-event population metrics
%       basename.rippStates.mat  - per-bout rate/density table (evt_states)
%       basename.rippMaps.mat    - per-event LFP maps
%       basename.rippSpks.mat    - per-unit spike stats + PETH (light)
%       basename.rippSpkMaps.mat - 3D spike raster [unit x event x bin] (heavy)
%       basename.rippSpkLfp.mat  - spike-LFP phase coupling
%
%   FILES SAVED (when flgNS = true; independent of flgSave):
%       basename.rip.evt         - NeuroScope events (start / peak / stop marks)
%       The lone output tagged 'rip' rather than the pipeline's 'ripp':
%       NeuroScope loads an event file only when the id is exactly three
%       characters, and refuses the file otherwise. See evt2ns.
%
%   DEPENDENCIES:
%       ripp_sigLoad, ripp_sigPrep, ripp_times, ripp_params, spklfp_phase,
%       guiPath; and the shared event layer (lfp/events): evt_files,
%       evt_pickCh, evt_spkPrep, evt_boutTimes, evt_emgScore, evt_spkGain,
%       evt_qa, evt_subset, evt_ctrlTimes, evt_states, evt_maps, evt_spks,
%       evt_saveSpks, evt_plotSpks, evt2ns.
%
%   HISTORY:
%       Updated: 260706 (flgForce-only re-run mirroring ed_wrapper; shared spine
%                via evt_* helpers; populate .info; drop the steps modes and the
%                Neuroscope export).
%       Updated: 260715 (restore the Neuroscope export as flgNS, now a .evt file
%                via evt2ns, placed after detect-or-load so it can be refreshed
%                post-curation).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'thr', [1, 3.5, 2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 20, 10], @isnumeric);
addParameter(p, 'thrEmg', 2, @isnumeric);
addParameter(p, 'win', [0, Inf], @isnumeric);
addParameter(p, 'passband', [80 250], @isnumeric);
addParameter(p, 'detectMet', 3, @isnumeric);
addParameter(p, 'zMet', 'nrem', @ischar);
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgNS', false, @islogical);
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

basepath  = p.Results.basepath;
rippCh    = p.Results.rippCh;
thr       = p.Results.thr;
limDur    = p.Results.limDur;
thrEmg    = p.Results.thrEmg;
win       = p.Results.win;
passband  = p.Results.passband;
detectMet = p.Results.detectMet;
zMet      = p.Results.zMet;
mapDur    = p.Results.mapDur;
bit2uv    = p.Results.bit2uv;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgNS     = p.Results.flgNS;
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
files = evt_files(basepath, basename, 'ripp');
files.phase = fullfile(basepath, [basename, '.rippSpkLfp.mat']);

if verbose, fprintf('[RIPP]: Session %s\n', basename); end

% Session + data. Spikes / units / states are optional and guarded downstream.
v = basepaths2vars('basepaths', {basepath}, ...
    'vars', {'session', 'spikes', 'spktimes', 'sleep_states', 'units'});
fsSpk = v.session.extracellular.sr;
if isinf(win(2)), sigDur = Inf; else, sigDur = win(2) - win(1); end

% Spike times (window-relative SU + pooled MUA) and bout times (both cheap;
% empty when the data is absent, which makes the dependent analyses skip).
[spkTimes, muTimes, uType] = evt_spkPrep(v, win, sigDur, fsSpk);
hasSpikes = ~isempty(spkTimes);
[boutTimes, vldTimes, nremTimes] = evt_boutTimes(v, win, sigDur);

%% ========================================================================
%  DETECT OR LOAD
%  ========================================================================
if isfile(files.evt) && ~flgForce
    if verbose, fprintf('[RIPP]: Loading existing %s.ripp.mat\n', basename); end
    S = load(files.evt, 'ripp');
    ripp = S.ripp;

else
    % ---- Signal ---------------------------------------------------------
    if verbose, fprintf('[RIPP]: Loading + pre-processing signal...\n'); end
    if isempty(rippCh), rippCh = evt_pickCh(v.session); end
    [lfp, emg, fs] = ripp_sigLoad(basepath, 'win', win, ...
        'session', v.session, 'basename', basename, ...
        'rippCh', rippCh, 'bit2uv', bit2uv);
    rippSig = ripp_sigPrep(lfp, fs, 'detectMet', detectMet, ...
        'passband', passband, 'zMet', zMet, 'nremTimes', nremTimes);

    % ---- Detect + QA ----------------------------------------------------
    if verbose, fprintf('[RIPP]: Thresholding candidate events...\n'); end
    ripp = ripp_times(rippSig, fs, 'thr', thr, 'limDur', limDur);

    % QA metrics: EMG (event EMG vs the NREM baseline) + MUA spike-gain
    ripp.emg = evt_emgScore(emg, ripp.times, fs, 'baselineTimes', nremTimes);
    ripp.spkGain = evt_spkGain(muTimes, ripp.times);

    % QA filter: valid vigilance states + low EMG + positive gain. Absent data
    % relaxes a criterion (empty inTimes or a NaN metric -> pass).
    idxGood = evt_qa(ripp.peakTime, ...
        'inTimes', vldTimes, ...
        'metrics', [ripp.emg, ripp.spkGain], ...
        'ranges', {[-Inf, thrEmg], [0, Inf]}, ...
        'names', {'emg', 'gain'});
    if verbose
        fprintf('[RIPP]: QA kept %d / %d events.\n', ...
            sum(idxGood), numel(idxGood));
    end
    ripp = evt_subset(ripp, idxGood);

    % ---- Characterise ---------------------------------------------------
    ripp.ctrlTimes = evt_ctrlTimes(ripp.times, ...
        'vldTimes', vldTimes, 'flgPlot', false);
    [ripp.state, ~] = evt_states(ripp.times, ripp.peakTime, boutTimes, ...
        'basepath', basepath, 'flgPlot', flgPlot, 'flgSave', flgSave, ...
        'name', 'ripp', 'lbl', 'Ripple');
    ripp = ripp_params(rippSig, ripp);
    if verbose, fprintf('[RIPP]: Generating LFP maps...\n'); end
    rippMaps = evt_maps(rippSig, ripp.peakTime, fs, ...
        'mapDur', mapDur, 'flgSave', false);

    % Seed curation acceptance (guiPath reads .accepted)
    ripp.accepted = true(size(ripp.times, 1), 1);

    % ---- Spiking (SU/MU modulation + PETH; phase) -----------------------
    hasSpks = hasSpikes && ~isempty(ripp.times);
    if hasSpks
        if verbose, fprintf('[RIPP]: Analysing spikes...\n'); end
        rippSpks = evt_spks(spkTimes, muTimes, ripp.times, ...
            ripp.ctrlTimes, ripp.peakTime, 'unitType', uType, 'mapDur', mapDur);
        ripp.spks = rippSpks.events;          % per-event metrics live on ripp
        rippSpks = rmfield(rippSpks, 'events');

        if verbose, fprintf('[RIPP]: Calculating spk-LFP phase...\n'); end
        spkLfp = spklfp_phase(rippSig.filt, spkTimes, fs, ...
            'lfpTimes', ripp.times, 'nPerms', 0);
    end

    % ---- Finalise: absolute time + provenance ---------------------------
    ripp.times    = ripp.times + win(1);
    ripp.peakTime = ripp.peakTime + win(1);

    % Merge provenance onto .info (ripp_times already set .fs/.thr/.limDur,
    % which ripp_params reads - assign per-field so they survive).
    ripp.info.basename  = basename;
    ripp.info.rippCh    = rippCh;
    ripp.info.passband  = passband;
    ripp.info.detectMet = detectMet;
    ripp.info.zMet      = zMet;
    ripp.info.win       = win;
    ripp.info.runtime   = datetime('now');

    % ---- Save -----------------------------------------------------------
    if flgSave
        if verbose, fprintf('[RIPP]: Saving output files...\n'); end
        save(files.evt, 'ripp', '-v7.3');
        save(files.maps, 'rippMaps', '-v7.3');
        if hasSpks
            evt_saveSpks(rippSpks, basepath, basename, 'ripp');
            save(files.phase, 'spkLfp', '-v7.3');
        end
    end

    % ---- Plot -----------------------------------------------------------
    if flgPlot && hasSpks
        if verbose, fprintf('[RIPP]: Generating summary plot...\n'); end
        evt_plotSpks(rippSpks, 'basepath', basepath, 'flgSaveFig', true, ...
            'name', 'ripp', 'lbl', 'Ripple');
    end
end

%% ========================================================================
%  NEUROSCOPE
%  ========================================================================
% Both branches leave ripp in absolute time, so the export reads the same on a
% fresh detection (.accepted all-true) and on a re-run over a curated session
% (the GUI's mask) - rejected events stay in the file under their own palette
% entry. A .ripp.mat predating the canonical schema carries no .accepted; an
% empty mask there labels every event as accepted rather than erroring.
if flgNS
    if verbose, fprintf('[RIPP]: Writing NeuroScope events...\n'); end
    accepted = [];
    if isfield(ripp, 'accepted'), accepted = ripp.accepted; end

    % 'rip', not 'ripp': NeuroScope rejects an event id that is not exactly
    % three characters, so this one file cannot follow the pipeline's tag.
    evt2ns(ripp.times, ripp.peakTime, 'basepath', basepath, ...
        'basename', basename, 'fileTag', 'rip', 'lbl', 'Ripple', ...
        'accepted', accepted);
end

%% ========================================================================
%  CURATE
%  ========================================================================
if flgCurate
    if isfile(files.evt)
        if verbose, fprintf('[RIPP]: Launching curation GUI...\n'); end
        guiPath(basepath, 'preset', 'Ripples', 'basename', basename);
    elseif verbose
        fprintf('[RIPP]: No %s.ripp.mat; set flgSave=true to curate.\n', ...
            basename);
    end
end

if verbose, fprintf('[RIPP]: Done (%s).\n', basename); end

end     % EOF
