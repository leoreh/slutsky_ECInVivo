function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Master function to execute the complete SWR analysis pipeline.
%
%   ripp = RIPP_WRAPPER(varargin)
%
%   SUMMARY:
%       Coordinates detection, analysis, and visualization of Sharp-Wave
%       Ripples (SWR). The body reads top-to-bottom as the pipeline stages:
%       1.  Setup: load session data; prepare spike times (evt_spkPrep).
%       2.  Signal: load the ripple LFP channel; filter + envelope + z-score
%           (ripp_sigPrep).
%       3.  Detect: threshold candidates (ripp_times); QA filter on state /
%           EMG / spike-gain (evt_qa); Neuroscope export (ripp2ns).
%       4.  Characterize: matched controls (evt_ctrlTimes); vigilance state
%           (evt_states); per-event params (ripp_params); LFP maps (evt_maps).
%       5.  Spiking: SU/MU modulation stats + PETHs, one call (evt_spks).
%       6.  Phase: spike-LFP coupling (spklfp_phase).
%       7.  Save / plot / curate (gui_curate).
%
%   INPUTS:
%       varargin - Parameter/Value pairs:
%           'basepath'   - (Char) Base directory of the session (default: pwd).
%           'rippCh'     - (Num)  Zero-indexed channel ID for ripple detection.
%                                 If empty, tries to load from session tags.
%           'thr'        - (Vec)  Detection thresholds [start, peak, cont, max, min_cont].
%           'limDur'     - (Vec)  Duration limits [min, max, inter, min_cont_dur] (ms).
%           'thrEmg'     - (Num)  EMG z-score pass threshold for QA. Default: 2.
%           'win'        - (Vec)  Time window to analyze [start end] (s). Default: [0 Inf].
%           'passband'   - (Vec)  Filtering frequency band [min max] (Hz). Default: [80 250].
%           'detectMet'  - (Num)  Detection method ID (see ripp_sigPrep). Default: 3.
%           'zMet'       - (Char) Z-scoring method ('adaptive', 'nrem'). Default: 'nrem'.
%           'mapDur'     - (Vec)  Window for PETH/Maps [pre post] (s). Default: [-0.1 0.1].
%           'bit2uv'     - (Num)  Conversion factor. Auto-detects TDT/Intan if empty.
%           'flgPlot'    - (Log)  Generate summary plots + viewers? Default: true.
%           'flgSave'    - (Log)  Save output .mat files? Default: false.
%           'flgCurate'  - (Log)  Launch the curation GUI (gui_curate)? Default: false.
%           'verbose'    - (Log)  Print progress steps? Default: true.
%           'steps'      - (Char) Execution mode: 'all', 'spks', 'phase'.
%                                 'all': Run full pipeline (default).
%                                 'spks': Skip signal/detect. Calc spikes/PETH.
%                                 'phase': Skip detect/spks. Calc phase.
%
%   OUTPUTS:
%       ripp         - (Struct) Combined structure of all analysis results:
%           .times       - Event start/end times.
%           .peakTime    - Event peak times.
%           .amp,.freq   - Event parameters.
%           .state       - Vigilance state per event.
%           .spkGain     - MUA gain per event.
%           .spks        - Per-event population spike metrics (frac/asym/com).
%
%   files saved (if flgSave=true):
%       basename.ripp.mat         - events + per-event population metrics
%       basename.rippStates.mat   - per-bout rate/density table
%       basename.rippMaps.mat     - per-event LFP maps (always producible)
%       basename.rippSpks.mat     - per-unit spike stats + per-unit PETH (light)
%       basename.rippSpkMaps.mat  - 3D spike raster [unit x event x bin] (heavy)
%       basename.rippSpkLfp.mat   - spike-LFP phase coupling
%
%   DEPENDENCIES:
%       ripp_sigPrep, ripp_times, ripp_params, ripp2ns, spklfp_phase, gui_curate;
%       and the shared event layer (lfp/events): evt_spkPrep, evt_emgScore,
%       evt_spkGain, evt_qa, evt_states, evt_ctrlTimes, evt_maps, evt_spks,
%       evt_plotSpks.
%
%   HISTORY:
%       Updated: 05 Jul 2026 (linear structure; algorithms in files; split
%                spike outputs into rippSpks + rippSpkMaps).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'thr', [1, 3.5, 2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 20, 10], @isnumeric);
addParameter(p, 'thrEmg', 2, @isnumeric);
addParameter(p, 'win', [0, Inf], @isnumeric);
addParameter(p, 'passband', [80 250], @isnumeric);
addParameter(p, 'detectMet', 3, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgCurate', false, @islogical);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'zMet', 'nrem', @ischar);
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'steps', 'all', @ischar);
parse(p, varargin{:});

basepath  = p.Results.basepath;
rippCh    = p.Results.rippCh;
thr       = p.Results.thr;
limDur    = p.Results.limDur;
thrEmg    = p.Results.thrEmg;
win       = p.Results.win;
passband  = p.Results.passband;
detectMet = p.Results.detectMet;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgCurate = p.Results.flgCurate;
bit2uv    = p.Results.bit2uv;
zMet      = p.Results.zMet;
mapDur    = p.Results.mapDur;
verbose   = p.Results.verbose;
steps     = p.Results.steps;

% Execution mode
doLoad   = ismember(steps, {'all', 'phase'});
doDetect = strcmpi(steps, 'all');
doSpks   = ismember(steps, {'all', 'spks'});
doPhase  = ismember(steps, {'all', 'phase'});

%% ========================================================================
%  SETUP
%  ========================================================================
[~, basename] = fileparts(basepath);
if verbose, fprintf('[RIPP]: Starting pipeline for %s...\n', basename); end

% Load session + data
vars = {'session', 'spikes', 'spktimes', 'sleep_states', 'units'};
if ~doDetect
    vars = [vars, {'ripp', 'rippMaps'}];
end
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

fs    = v.session.extracellular.srLfp;
fsSpk = v.session.extracellular.sr;
if isinf(win(2)), sigDur = Inf; else, sigDur = win(2) - win(1); end

% Spike preparation (shared): window-relative single-unit times, one pooled
% MUA vector (for QA spike-gain), and unit types. Absent data yields empty
% outputs, so spike analyses skip and the cell-type split falls back to "Global".
[spkTimes, muTimes, uType] = evt_spkPrep(v, win, sigDur, fsSpk);
hasSpikes = ~isempty(spkTimes);
if ~hasSpikes && verbose
    fprintf('[RIPP]: No sorted spikes; skipping spike analyses.\n');
end
doSpks  = doSpks  && hasSpikes;
doPhase = doPhase && hasSpikes;

% Sleep-state bout times (window-relative). NREM (4th cell) seeds the EMG
% baseline; valid states seed control matching. Missing states relax both.
try
    boutTimes = v.ss.bouts.times;
    boutTimes = cellfun(@(x) x - win(1), boutTimes, 'UniformOutput', false);
    boutTimes = cellfun(@(x) x(x(:,2)>0 & x(:,1)<sigDur, :), boutTimes, 'UniformOutput', false);
    nremTimes = boutTimes{4};
    vldTimes = vertcat(boutTimes{2}, boutTimes{3}, boutTimes{4});
catch
    warning('Could not load sleep states.');
    boutTimes = [];
    nremTimes = [];
    vldTimes = [];
end

% Output files
files.ripp    = fullfile(basepath, [basename, '.ripp.mat']);
files.maps    = fullfile(basepath, [basename, '.rippMaps.mat']);
files.spks    = fullfile(basepath, [basename, '.rippSpks.mat']);
files.spkMaps = fullfile(basepath, [basename, '.rippSpkMaps.mat']);
files.phase   = fullfile(basepath, [basename, '.rippSpkLfp.mat']);

%% ========================================================================
%  SIGNAL  (load the ripple LFP channel + filter / envelope / z-score)
%  ========================================================================
rippSig = [];
if doLoad
    if verbose, fprintf('[RIPP]: Loading + pre-processing signal...\n'); end

    nchans = v.session.extracellular.nChannels;

    % Ripple channel (from session tags unless given)
    if isempty(rippCh)
        if isfield(v.session.channelTags, 'Ripple')
            rippCh = v.session.channelTags.Ripple;
        else
            rippCh = 1;
        end
    end

    % Voltage conversion (auto: TDT vs Intan)
    if isempty(bit2uv)
        if round(v.session.extracellular.sr) == 24414
            bit2uv = 1;      % TDT / Tucker
        else
            bit2uv = 0.195;  % Intan
        end
    end

    % Load the channel(s), average if more than one
    fname = fullfile(basepath, [basename, '.lfp']);
    lfp = double(binary_load(fname, 'duration', sigDur, 'fs', fs, 'nCh', nchans, ...
        'start', win(1), 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv));
    if size(lfp, 2) > 1, lfp = mean(lfp, 2); end

    % Filter + envelope + z-score
    rippSig = ripp_sigPrep(lfp, fs, ...
        'detectMet', detectMet, 'passband', passband, ...
        'zMet', zMet, 'nremTimes', nremTimes);
end

%% ========================================================================
%  DETECT  (threshold -> QA filter -> Neuroscope)
%  ========================================================================
if doDetect

    % EMG (optional) for the QA criterion; missing / misfit EMG relaxes it.
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
        if length(emg) ~= length(rippSig.lfp)
            warning('ripp_wrapper:emgFit', ...
                'EMG (%d) and LFP (%d) lengths differ; skipping EMG QA.', ...
                length(emg), length(rippSig.lfp));
            emg = [];
        end
    end

    % Candidate events
    if verbose, fprintf('[RIPP]: Thresholding candidate events...\n'); end
    ripp = ripp_times(rippSig, fs, 'thr', thr, 'limDur', limDur);

    % QA metrics: EMG (event EMG vs the NREM baseline) + MUA spike-gain.
    ripp.emg = evt_emgScore(emg, ripp.times, fs, 'baselineTimes', nremTimes);
    ripp.spkGain = evt_spkGain(muTimes, ripp.times);

    % QA filter: keep events in valid vigilance states (inTimes), with low EMG
    % and positive spike-gain. Each criterion is skipped when its data is absent
    % (empty inTimes or a NaN metric -> pass).
    idxGood = evt_qa(ripp.peakTime, ...
        'inTimes', vldTimes, ...
        'metrics', [ripp.emg, ripp.spkGain], ...
        'ranges', {[-Inf, thrEmg], [0, Inf]}, ...
        'names', {'emg', 'gain'});
    if verbose
        fprintf('[RIPP]: QA kept %d / %d events.\n', sum(idxGood), numel(idxGood));
    end

    % Filter every per-event field by the QA mask (drop the failures)
    fnames = fieldnames(ripp);
    for iField = 1:numel(fnames)
        fn = fnames{iField};
        if size(ripp.(fn), 1) == numel(idxGood)
            ripp.(fn) = ripp.(fn)(idxGood, :);
        end
    end

    % Neuroscope export (absolute samples)
    if flgSave
        if verbose, fprintf('[RIPP]: Saving Neuroscope events...\n'); end
        rippSamps = round((ripp.times + win(1)) * fsSpk);
        peakSamps = round((ripp.peakTime + win(1)) * fsSpk);
        ripp2ns(rippSamps, peakSamps, 'basepath', basepath);
    end

    % Matched control intervals
    ripp.ctrlTimes = evt_ctrlTimes(ripp.times, 'vldTimes', vldTimes, 'flgPlot', false);

    % Vigilance state (final; saves rippStates)
    [ripp.state, ~] = evt_states(ripp.times, ripp.peakTime, boutTimes, ...
        'basepath', basepath, 'flgPlot', flgPlot, 'flgSave', true, ...
        'name', 'ripp', 'lbl', 'Ripple');

    % Per-event params + LFP maps
    ripp = ripp_params(rippSig, ripp);
    if verbose, fprintf('[RIPP]: Generating LFP maps...\n'); end
    rippMaps = evt_maps(rippSig, ripp.peakTime, fs, 'mapDur', mapDur, 'flgSave', false);

    % Seed curation acceptance (gui_curate reads .accepted)
    ripp.accepted = true(size(ripp.times, 1), 1);

else
    if verbose, fprintf('[RIPP]: Skipping detection (loading from file)...\n'); end
    if isfield(v, 'ripp'), ripp = v.ripp; end
    if isfield(v, 'rippMaps'), rippMaps = v.rippMaps; end
end

%% ========================================================================
%  SPIKE STATS  (per-unit stats + population metrics + PETH, one call)
%  ========================================================================
if doSpks
    if verbose, fprintf('[RIPP]: Analyzing spikes...\n'); end
    rippSpks = evt_spks(spkTimes, muTimes, ...
        ripp.times, ripp.ctrlTimes, ripp.peakTime, ...
        'unitType', uType, 'mapDur', mapDur);

    % Per-event population metrics live on the ripp struct
    ripp.spks = rippSpks.events;
    rippSpks = rmfield(rippSpks, 'events');
end

%% ========================================================================
%  PHASE COUPLING
%  ========================================================================
if doPhase
    if verbose, fprintf('[RIPP]: Calculating Spk-LFP phase...\n'); end
    spkLfp = spklfp_phase(rippSig.filt, spkTimes, fs, ...
        'lfpTimes', ripp.times, 'nPerms', 0);
end

%% ========================================================================
%  SAVE
%  ========================================================================
% Absolute times for saving / Neuroscope
ripp.times = ripp.times + win(1);
ripp.peakTime = ripp.peakTime + win(1);

if flgSave
    if verbose, fprintf('[RIPP]: Saving output files...\n'); end

    if doDetect || doSpks
        save(files.ripp, 'ripp', '-v7.3');
        save(files.maps, 'rippMaps', '-v7.3');
    end

    % Split the spike outputs: heavy 3D raster in its own file, light per-unit
    % stats + PETH in rippSpks (so the manuscript loader stays fast).
    if doSpks
        rippSpkMaps = rippSpks.maps;            % .su/.mu, each .evt/.ctrl
        rippSpkMaps.tstamps = rippSpks.tstamps;
        save(files.spkMaps, 'rippSpkMaps', '-v7.3');

        % Write rippSpks without the raster; keep the in-memory copy for plots.
        spksOut = struct('rippSpks', rmfield(rippSpks, 'maps'));
        save(files.spks, '-struct', 'spksOut', '-v7.3');
    end

    if doPhase
        save(files.phase, 'spkLfp', '-v7.3');
    end
end

%% ========================================================================
%  PLOT
%  ========================================================================
if flgPlot && doSpks
    if verbose, fprintf('[RIPP]: Generating summary plot...\n'); end
    evt_plotSpks(rippSpks, 'basepath', basepath, 'flgSaveFig', true, ...
        'name', 'ripp', 'lbl', 'Ripple');
end

%% ========================================================================
%  CURATE
%  ========================================================================
if flgCurate
    if isfile(files.ripp)
        if verbose, fprintf('[RIPP]: Launching curation GUI...\n'); end
        gui_curate(basepath, 'preset', 'Ripples');
    elseif verbose
        fprintf('[RIPP]: No %s on disk; set flgSave=true to curate.\n', [basename, '.ripp.mat']);
    end
end

if verbose, fprintf('[RIPP]: Pipeline completed for %s.\n', basename); end

end     % EOF
